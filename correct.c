#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include "bfc.h"

/**************************
 * Sequence struct for ec *
 **************************/

#include "kvec.h"

typedef struct { // NOTE: unaligned memory
	uint8_t b:3, q:1, ob:3, oq:1;
	uint8_t dummy;
	uint16_t lcov:6, hcov:6, solid_end:1, high_end:1, ec:1, absent:1;
	int i;
} ecbase_t;

typedef kvec_t(ecbase_t) ecseq_t;

static int bfc_seq_conv(const char *s, const char *q, int qthres, ecseq_t *seq, int b_from_q)
{
	int i, l;
	l = strlen(s);
	kv_resize(ecbase_t, *seq, l);
	seq->n = l;
	for (i = 0; i < l; ++i) {
		ecbase_t *c = &seq->a[i];
		c->b = c->ob = b_from_q && q && q[i] - 33 <= 5? q[i] - 34 : seq_nt6_table[(int)s[i]] - 1;
		c->q = c->oq = !q? 1 : q[i] - 33 >= qthres? 1 : 0;
		if (c->b > 3) c->q = c->oq = 0;
		c->i = i;
	}
	return l;
}

static inline ecbase_t ecbase_comp(const ecbase_t *b)
{
	ecbase_t r = *b;
	r.b = b->b < 4? 3 - b->b : 4;
	r.ob = b->ob < 4? 3 - b->ob : 4;
	return r;
}

static void bfc_seq_revcomp(ecseq_t *seq)
{
	int i;
	for (i = 0; i < seq->n>>1; ++i) {
		ecbase_t tmp;
		tmp = ecbase_comp(&seq->a[i]);
		seq->a[i] = ecbase_comp(&seq->a[seq->n - 1 - i]);
		seq->a[seq->n - 1 - i] = tmp;
	}
	if (seq->n&1) seq->a[i] = ecbase_comp(&seq->a[i]);
}

/***************************
 * Independent ec routines *
 ***************************/

int bfc_ec_greedy_k(int k, int mode, const bfc_kmer_t *x, const bfc_ch_t *ch)
{
	int i, j, max = 0, max_ec = -1, max2 = 0;
	for (i = 0; i < k; ++i) {
		int c = (x->x[1]>>i&1)<<1 | (x->x[0]>>i&1);
		for (j = 0; j < 4; ++j) {
			bfc_kmer_t y = *x;
			int ret;
			if (j == c) continue;
			bfc_kmer_change(k, y.x, i, j);
			ret = bfc_ch_kmer_occ(ch, &y);
			if (ret < 0) continue;
			if ((max&0xff) < (ret&0xff)) max2 = max, max = ret, max_ec = i<<2 | j;
			else if ((max2&0xff) < (ret&0xff)) max2 = ret;
		}
	}
	return (max&0xff) * 3 > mode && (max2&0xff) < 3? max_ec : -1;
}

int bfc_ec_first_kmer(int k, const ecseq_t *s, int start, bfc_kmer_t *x)
{
	int i, l;
	*x = bfc_kmer_null;
	for (i = start, l = 0; i < s->n; ++i) {
		ecbase_t *c = &s->a[i];
		if (c->b < 4) {
			bfc_kmer_append(k, x->x, c->b);
			if (++l == k) break;
		} else l = 0, *x = bfc_kmer_null;
	}
	return i;
}

void bfc_ec_kcov(int k, int min_occ, ecseq_t *s, const bfc_ch_t *ch)
{
	int i, l, r, j;
	bfc_kmer_t x = bfc_kmer_null;
	for (i = l = 0; i < s->n; ++i) {
		ecbase_t *c = &s->a[i];
		c->high_end = c->solid_end = c->lcov = c->hcov = 0;
		if (c->b < 4) {
			bfc_kmer_append(k, x.x, c->b);
			if (++l >= k) {
				if ((r = bfc_ch_kmer_occ(ch, &x)) >= 0) {
					if ((r>>8&0x3f) >= min_occ+1) c->high_end = 1;
					if ((r&0xff) >= min_occ) {
						c->solid_end = 1;
						for (j = i - k + 1; j <= i; ++j)
							++s->a[j].lcov, s->a[j].hcov += c->high_end;
					}
				}
			}
		} else l = 0, x = bfc_kmer_null;
	}
}

uint64_t bfc_ec_best_island(int k, const ecseq_t *s)
{ // IMPORTANT: call bfc_ec_kcov() before calling this function!
	int i, l, max, max_i;
	for (i = k - 1, max = l = 0, max_i = -1; i < s->n; ++i) {
		if (!s->a[i].solid_end) {
			if (l > max) max = l, max_i = i;
			l = 0;
		} else ++l;
	}
	if (l > max) max = l, max_i = i;
	return max > 0? (uint64_t)(max_i - max - k + 1) << 32 | max_i : 0;
}

/********************
 * Correct one read *
 ********************/

#include "ksort.h"

#define ECCODE_MISC      1
#define ECCODE_MANY_N    2
#define ECCODE_NO_SOLID  3
#define ECCODE_UNCORR_N  4
#define ECCODE_MANY_FAIL 5

typedef struct {
	uint32_t ec_code:3, brute:1, n_ec:14, n_ec_high:14;
	uint32_t n_absent:22, rf_code:2, max_heap:8;
} ecstat_t;

typedef struct {
	uint8_t ec:1, ec_high:1, absent:1, absent_high:1, b:4;
} bfc_penalty_t;

typedef struct {
	int tot_pen;
	int i; // base position
	int k; // position in the stack
	int32_t ecpos_high[BFC_EC_HIST_HIGH];
	int32_t ecpos[BFC_EC_HIST];
	bfc_kmer_t x;
} echeap1_t;

typedef struct {
	int parent, i, tot_pen;
	uint8_t b;
	bfc_penalty_t pen;
	uint16_t cnt;
} ecstack1_t;

typedef struct {
	const bfc_opt_t *opt;
	const bfc_ch_t *ch;
	kvec_t(echeap1_t) heap;
	kvec_t(ecstack1_t) stack;
	ecseq_t seq, tmp, ec[2];
	int mode;
	ecstat_t ori_st;
} bfc_ec1buf_t;

#define heap_lt(a, b) ((a).tot_pen > (b).tot_pen)
KSORT_INIT(ec, echeap1_t, heap_lt)

static bfc_ec1buf_t *ec1buf_init(const bfc_opt_t *opt, const bfc_ch_t *ch)
{
	bfc_ec1buf_t *e;
	e = calloc(1, sizeof(bfc_ec1buf_t));
	e->opt = opt, e->ch = ch;
	return e;
}

static void ec1buf_destroy(bfc_ec1buf_t *e)
{	
	free(e->heap.a); free(e->stack.a); free(e->seq.a); free(e->tmp.a); free(e->ec[0].a); free(e->ec[1].a);
	free(e);
}

#define weighted_penalty(o, p) ((o)->w_ec * (p).ec + (o)->w_ec_high * (p).ec_high + (o)->w_absent * (p).absent + (o)->w_absent_high * (p).absent_high)

static void buf_update(bfc_ec1buf_t *e, const echeap1_t *prev, bfc_penalty_t pen, int cnt)
{
	ecstack1_t *q;
	echeap1_t *r;
	const bfc_opt_t *o = e->opt;
	int b = pen.b;
	// update stack
	kv_pushp(ecstack1_t, e->stack, &q);
	q->parent = prev->k;
	q->i = prev->i;
	q->b = b;
	q->pen = pen;
	q->cnt = cnt > 0? cnt&0xff : 0;
	q->tot_pen = prev->tot_pen + weighted_penalty(o, pen);
	// update heap
	kv_pushp(echeap1_t, e->heap, &r);
	r->i = prev->i + 1;
	r->k = e->stack.n - 1;
	r->x = prev->x;
	if (pen.ec_high) {
		memcpy(r->ecpos_high + 1, prev->ecpos_high, (BFC_EC_HIST_HIGH - 1) * 4);
		r->ecpos_high[0] = prev->i;
	} else memcpy(r->ecpos_high, prev->ecpos_high, BFC_EC_HIST_HIGH * 4);
	if (pen.ec) {
		memcpy(r->ecpos + 1, prev->ecpos, (BFC_EC_HIST - 1) * 4);
		r->ecpos[0] = prev->i;
	} else memcpy(r->ecpos, prev->ecpos, BFC_EC_HIST * 4);
	r->tot_pen = q->tot_pen;
	bfc_kmer_append(e->opt->k, r->x.x, b);
	if (bfc_verbose >= 4)
		fprintf(stderr, "     <= base:%c penalty:%d\n", pen.ec? "acgtn"[b] : "ACGTN"[b], r->tot_pen);
	ks_heapup_ec(e->heap.n, e->heap.a);
}

static int buf_backtrack(ecstack1_t *s, int end, const ecseq_t *seq, ecseq_t *path)
{
	int i, n_absent = 0;
	kv_resize(ecbase_t, *path, seq->n);
	path->n = seq->n;
	while (end >= 0) {
		if ((i = s[end].i) < seq->n) {
			path->a[i].b = s[end].b;
			path->a[i].ec = s[end].pen.ec;
			path->a[i].absent = s[end].pen.absent;
			n_absent += s[end].pen.absent;
		}
		end = s[end].parent;
	}
	return n_absent;
}

static int bfc_ec1dir(bfc_ec1buf_t *e, const ecseq_t *seq, ecseq_t *ec, int start, int end, int *max_heap)
{
	echeap1_t z;
	int i, l, rv = -1, path[BFC_MAX_PATHS], n_paths = 0, min_path = -1, min_path_pen = INT_MAX, n_failures = 0;
	assert(end <= seq->n && end - start >= e->opt->k);
	if (bfc_verbose >= 4) fprintf(stderr, "* bfc_ec1dir(): len:%ld start:%d end:%d\n", seq->n, start, end);
	e->heap.n = e->stack.n = 0;
	*max_heap = 0;
	memset(&z, 0, sizeof(echeap1_t));
	kv_resize(ecbase_t, *ec, seq->n);
	ec->n = seq->n;
	for (z.i = start, l = 0; z.i < end; ++z.i) {
		int c = seq->a[z.i].b;
		if (c < 4) {
			if (++l == e->opt->k) break;
			bfc_kmer_append(e->opt->k, z.x.x, c);
		} else l = 0, z.x = bfc_kmer_null;
	}
	assert(z.i < end); // before calling this function, there must be at least one solid k-mer
	z.k = -1;
	for (i = 0; i < BFC_EC_HIST; ++i) z.ecpos[i] = -1;
	for (i = 0; i < BFC_EC_HIST_HIGH; ++i) z.ecpos_high[i] = -1;
	kv_push(echeap1_t, e->heap, z);
	for (i = 0; i < seq->n; ++i) ec->a[i].b = seq->a[i].b, ec->a[i].ob = seq->a[i].ob;
	// exhaustive error correction
	while (1) {
		int stop = 0;
		*max_heap = *max_heap > 255? 255 : *max_heap > e->heap.n? *max_heap : e->heap.n;
		if (e->heap.n == 0) { // may happen when there is an uncorrectable "N"
			rv = -2;
			break;
		}
		z = e->heap.a[0];
		e->heap.a[0] = kv_pop(e->heap);
		ks_heapdown_ec(0, e->heap.n, e->heap.a);
		if (bfc_verbose >= 4)
			fprintf(stderr, "  => pos:%d stack_size:%ld heap_size:%ld penalty:%d last_base:%c ecpos_high:[%d,%d] ecpos:[%d,%d,%d,%d,%d]\n",
					z.i, e->stack.n, e->heap.n, z.tot_pen, "ACGT"[(z.x.x[1]&1)<<1|(z.x.x[0]&1)], z.ecpos_high[0], z.ecpos_high[1],
					z.ecpos[0], z.ecpos[1], z.ecpos[2], z.ecpos[3], z.ecpos[4]);
		if (min_path >= 0 && z.tot_pen > min_path_pen + e->opt->max_path_diff) break;
		if (z.i - end > e->opt->max_end_ext) stop = 1;
		if (!stop) {
			ecbase_t *c = z.i < seq->n? &seq->a[z.i] : 0;
			int b, os = -1, fixed = 0, other_ext = 0, n_added = 0, added_cnt[4];
			bfc_penalty_t added[4];
			// test if the read extension alone is enough
			if (z.i > end) fixed = 1;
			if (c && c->b < 4) { // A, C, G or T
				bfc_kmer_t x = z.x;
				bfc_kmer_append(e->opt->k, x.x, c->b);
				os = bfc_ch_kmer_occ(e->ch, &x);
				if (c->q && (os&0xff) >= e->opt->min_cov + 1 && c->lcov >= e->opt->min_cov + 1) fixed = 1;
				else if (c->hcov > e->opt->k * .75) fixed = 1;
				if (bfc_verbose >= 4) {
					fprintf(stderr, "     Original base:%c qual:%d fixed:%d count:", "ACGTN"[c->b], c->q, fixed);
					if (os >= 0) fprintf(stderr, "%d,%d\n", os&0xff, os>>8&0x3f);
					else fprintf(stderr, "-1,-1\n");
				}
			}
			// extension
			for (b = 0; b < 4; ++b) {
				bfc_penalty_t pen;
				if (fixed && c && b != c->b) continue;
				if (c == 0 || b != c->b) {
					int s;
					bfc_kmer_t x = z.x;
					if (c) { // not over the end
						if (c->q && z.ecpos_high[BFC_EC_HIST_HIGH-1] >= 0 && z.i - z.ecpos_high[BFC_EC_HIST_HIGH-1] < e->opt->win_multi_ec) continue; // no close highQ corrections
						if (z.ecpos[BFC_EC_HIST-1] >= 0 && z.i - z.ecpos[BFC_EC_HIST-1] < e->opt->win_multi_ec) continue; // no clustered corrections
					}
					bfc_kmer_append(e->opt->k, x.x, b);
					s = bfc_ch_kmer_occ(e->ch, &x);
					if (bfc_verbose >= 4 && s >= 0)
						fprintf(stderr, "     Alternative k-mer count: %c,%d:%d\n", "ACGTN"[b], s&0xff, s>>8&0x3f);
					if (s < 0 || (s&0xff) < e->opt->min_cov) continue; // not solid
					//if (os >= 0 && (s&0xff) - (os&0xff) < 2) continue; // not sufficiently better than the read path
					pen.ec = c && c->b < 4? 1 : 0;
					pen.ec_high = pen.ec? c->oq : 0;
					pen.absent = 0;
					pen.absent_high = ((s>>8&0xff) < e->opt->min_cov);
					pen.b = b;
					added_cnt[n_added] = s;
					added[n_added++] = pen;
					++other_ext;
				} else {
					pen.ec = pen.ec_high = 0;
					pen.absent = (os < 0 || (os&0xff) < e->opt->min_cov);
					pen.absent_high = (os < 0 || (os>>8&0xff) < e->opt->min_cov);
					pen.b = b;
					added_cnt[n_added] = os;
					added[n_added++] = pen;
				}
			} // ~for(b)
			if (fixed == 0 && other_ext == 0) ++n_failures;
			if (n_failures > seq->n * 2) {
				if (bfc_verbose >= 4) fprintf(stderr, "  !! too many unsuccessful attempts\n");
				rv = -3;
				break;
			}
			if (c || n_added == 1) {
				if (n_added > 1 && e->heap.n > e->opt->max_heap) { // to prevent heap explosion
					int min_b = -1, min = INT_MAX;
					for (b = 0; b < n_added; ++b) {
						int t = weighted_penalty(e->opt, added[b]);
						if (min > t) min = t, min_b = b;
					}
					buf_update(e, &z, added[min_b], added_cnt[min_b]);
				} else {
					for (b = 0; b < n_added; ++b)
						buf_update(e, &z, added[b], added_cnt[b]);
				}
			} else {
				if (n_added == 0)
					e->stack.a[z.k].tot_pen += e->opt->w_absent * (e->opt->max_end_ext - (z.i - end));
				stop = 1;
			}
		} // ~if(!stop)
		if (stop) {
			if (e->stack.a[z.k].tot_pen < min_path_pen)
				min_path_pen = e->stack.a[z.k].tot_pen, min_path = n_paths;
			path[n_paths++] = z.k;
			if (bfc_verbose >= 4) fprintf(stderr, "  @@ n_paths=%d penalty=%d\n", n_paths, e->stack.a[z.k].tot_pen);
			if (n_paths == BFC_MAX_PATHS) break;
		}
	} // ~while(1)
	// backtrack
	if (n_paths == 0) return rv;
	assert(min_path >= 0 && min_path < n_paths && e->stack.a[path[min_path]].tot_pen == min_path_pen);
	rv = buf_backtrack(e->stack.a, path[min_path], seq, ec);
	for (i = 0; i < ec->n; ++i) // mask out uncorrected regions
		if (i < start + e->opt->k || i >= end) ec->a[i].b = 4;
	if (bfc_verbose >= 4) {
		fprintf(stderr, "* %d path(s); lowest penalty: %d\n  ", n_paths, min_path_pen);
		for (i = 0; i < ec->n; ++i) fputc((seq->a[i].b == ec->a[i].b? "ACGTN":"acgtn")[ec->a[i].b], stderr);
		fputc('\n', stderr);
	}
	return rv;
}

ecstat_t bfc_ec1(bfc_ec1buf_t *e, char *seq, char *qual)
{
	int i, start = 0, end = 0, n_n = 0, rv[2], max_heap[2];
	uint64_t r;
	ecstat_t s;

	s.ec_code = ECCODE_MISC, s.brute = 0, s.n_ec = s.n_ec_high = 0, s.n_absent = s.max_heap = 0;
	s.rf_code = e->opt->refine_ec? 1 : 0;
	bfc_seq_conv(seq, qual, e->opt->q, &e->seq, e->opt->refine_ec);
	for (i = 0; i < e->seq.n; ++i)
		if (e->seq.a[i].ob > 3) ++n_n;
	if (n_n > e->seq.n * .05) {
		s.ec_code = ECCODE_MANY_N;
		return s;
	}
	bfc_ec_kcov(e->opt->k, e->opt->min_cov, &e->seq, e->ch);
	r = bfc_ec_best_island(e->opt->k, &e->seq);
	if (r == 0) { // no solid k-mer
		bfc_kmer_t x;
		int ec = -1;
		while ((end = bfc_ec_first_kmer(e->opt->k, &e->seq, start, &x)) < e->seq.n) {
			ec = bfc_ec_greedy_k(e->opt->k, e->mode, &x, e->ch);
			if (ec >= 0) break;
			if (end + (e->opt->k>>1) >= e->seq.n) break;
			start = end - (e->opt->k>>1);
		}
		if (ec >= 0) {
			e->seq.a[end - (ec>>2)].b = ec&3;
			++end; start = end - e->opt->k;
			s.brute = 1;
		} else {
			s.ec_code = ECCODE_NO_SOLID;
			return s;
		}
	} else start = r>>32, end = (uint32_t)r;
	if (bfc_verbose >= 4)
		fprintf(stderr, "* Longest solid island: [%d,%d)\n", start, end);
	if ((rv[0] = bfc_ec1dir(e, &e->seq, &e->ec[0], start, e->seq.n, &max_heap[0])) < 0) {
		s.ec_code = rv[0] == -2? ECCODE_UNCORR_N : rv[0] == -3? ECCODE_MANY_FAIL : ECCODE_MISC;
		return s;
	}
	bfc_seq_revcomp(&e->seq);
	if ((rv[1] = bfc_ec1dir(e, &e->seq, &e->ec[1], e->seq.n - end, e->seq.n, &max_heap[1])) < 0) {
		s.ec_code = rv[1] == -2? ECCODE_UNCORR_N : rv[1] == -3? ECCODE_MANY_FAIL : ECCODE_MISC;
		return s;
	}
	s.max_heap = max_heap[0] > max_heap[1]? max_heap[0] : max_heap[1];
	s.ec_code = 0, s.n_absent = rv[0] + rv[1];
	bfc_seq_revcomp(&e->ec[1]);
	bfc_seq_revcomp(&e->seq);
	if (e->opt->refine_ec && e->ori_st.ec_code == 0 && s.n_absent > e->ori_st.n_absent) {
		s = e->ori_st;
		s.rf_code = 2;
		return s;
	}
	for (i = 0; i < e->seq.n; ++i) {
		ecbase_t *c = &e->seq.a[i];
		if (e->ec[0].a[i].b == e->ec[1].a[i].b)
			c->b = e->ec[0].a[i].b > 3? e->seq.a[i].b : e->ec[0].a[i].b;
		else if (e->ec[1].a[i].b > 3) c->b = e->ec[0].a[i].b;
		else if (e->ec[0].a[i].b > 3) c->b = e->ec[1].a[i].b;
		else c->b = e->seq.a[i].ob;
	}
	for (i = 0; i < e->seq.n; ++i) {
		int is_diff = !(e->seq.a[i].b == e->seq.a[i].ob);
		if (is_diff) {
			++s.n_ec;
			if (e->seq.a[i].q) ++s.n_ec_high;
		}
		seq[i] = (is_diff? "acgtn" : "ACGTN")[e->seq.a[i].b];
		if (qual) qual[i] = is_diff? 34 + e->seq.a[i].ob : "+?"[e->seq.a[i].q];
	}
	if (bfc_verbose >= 4) {
		bfc_ec_kcov(e->opt->k, e->opt->min_cov, &e->seq, e->ch);
		fprintf(stderr, "* ec_code:%d n_ec:%d n_ec_high:%d\n  ", s.ec_code, s.n_ec, s.n_ec_high);
		for (i = 0; i < e->seq.n; ++i)
			fputc((e->seq.a[i].b == e->seq.a[i].ob? "ACGTN" : "acgtn")[e->seq.a[i].b], stderr);
		fprintf(stderr, "\n  ");
		for (i = 0; i < e->seq.n; ++i)
			fputc('0' + (int)(10. * e->seq.a[i].lcov / e->opt->k + .499), stderr);
		fputc('\n', stderr);
	}
	if (e->opt->refine_ec) s.rf_code = 3;
	return s;
}

/******************
 * K-mer trimming *
 ******************/

static uint64_t max_streak(int k, const bfc_bf_t *bf, const bseq1_t *s)
{
	int i, l;
	uint64_t max = 0, t = 0;
	bfc_kmer_t x = bfc_kmer_null;
	for (i = l = 0; i < s->l_seq; ++i) {
		int c = seq_nt6_table[(uint8_t)s->seq[i]] - 1;
		if (c < 4) { // not an ambiguous base
			bfc_kmer_append(k, x.x, c);
			if (++l >= k) { // ok, we have a k-mer now
				uint64_t hash, y[2];
				hash = bfc_kmer_hash(k, x.x, y);
				if (bfc_bf_get(bf, hash) == bf->n_hashes) t += 1ULL<<32; // in bloom filter
				else t = i + 1;
			} else t = i + 1;
		} else l = 0, x = bfc_kmer_null, t = i + 1;
		max = max > t? max : t;
	}
	return max;
}

/********************
 * Error correction *
 ********************/

typedef struct {
	const bfc_opt_t *opt;
	bseq_file_t *ks;
	const bfc_bf_t *bf;
	bfc_ec1buf_t **e;
	int64_t n_processed;
} ec_shared_t;

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	ec_shared_t *es;
} ec_step_t;

static inline ecstat_t parse_stats(char *str)
{
	ecstat_t s;
	char *p = str;
	s.ec_code = strtol(p, &p, 10);
	s.rf_code = 1;
	if (s.ec_code == 0) {
		s.n_absent = strtol(p + 1, &p, 10);
		s.max_heap = strtol(p + 1, &p, 10);
		s.brute = strtol(p + 1, &p, 10);
		s.n_ec = strtol(p + 1, &p, 10);
		s.n_ec_high = strtol(p + 1, &p, 10);
	} else s.n_absent = s.max_heap = s.brute = s.n_ec = s.n_ec_high = 0;
	return s;
}

static void worker_ec(void *_data, long k, int tid)
{
	ec_step_t *data = (ec_step_t*)_data;
	ec_shared_t *es = data->es;
	bseq1_t *s = &data->seqs[k];
	if (!es->opt->filter_mode) {
		ecstat_t st;
		bfc_ec1buf_t *e = es->e[tid];
		if (bfc_verbose >= 4) fprintf(stderr, "* Processing read '%s'...\n", s->name);
		if (es->opt->refine_ec && s->comment && strncmp(s->comment, "ec:Z:", 5) == 0) {
			e->ori_st = parse_stats(s->comment + 5);
			if (e->ori_st.ec_code == 0 && e->ori_st.max_heap < 50)
				return;
		}
		if (s->comment) {
			free(s->comment);
			s->comment = 0;
		}
		st = bfc_ec1(e, s->seq, s->qual);
		s->aux = st.n_ec<<18 | st.n_ec_high<<4 | st.brute<<3 | st.ec_code;
		s->aux2 = st.n_absent << 10 | st.rf_code << 8 | st.max_heap;
	} else {
		uint64_t max;
		max = max_streak(es->opt->k, es->bf, s);
		if (max>>32 && (double)((max>>32) + es->opt->k) / s->l_seq > es->opt->min_frac) {
			int start = (uint32_t)max, end = start + (max>>32);
			start -= es->opt->k - 1;
			assert(start >= 0 && end <= s->l_seq);
			memmove(s->seq, s->seq + start, end - start);
			s->l_seq = end - start;
			s->seq[s->l_seq] = 0;
			if (s->qual) {
				memmove(s->qual, s->qual + start, s->l_seq);
				s->qual[s->l_seq] = 0;
			}
			s->aux = 0;
		} else s->aux = 1;
	}
}

void *bfc_ec_cb(void *shared, int step, void *_data)
{
	ec_shared_t *es = (ec_shared_t*)shared;
	if (step == 0) {
		ec_step_t *ret;
		int keep_comment = (es->opt->filter_mode || es->opt->refine_ec);
		ret = calloc(1, sizeof(ec_step_t));
		ret->seqs = bseq_read(es->ks, es->opt->chunk_size, keep_comment, &ret->n_seqs);
		ret->es = es;
		fprintf(stderr, "[M::%s] read %d sequences\n", __func__, ret->n_seqs);
		if (ret->seqs) return ret;
		else free(ret);
	} else if (step == 1) {
		ec_step_t *data = (ec_step_t*)_data;
		kt_for(es->opt->n_threads, worker_ec, data, data->n_seqs);
		fprintf(stderr, "[M::%s @%.1f*%.1f%%] processed %d sequences\n", __func__, realtime() - bfc_real_time,
				100.*cputime()/(realtime()-bfc_real_time+1e-6), data->n_seqs);
		return data;
	} else if (step == 2) {
		ec_step_t *data = (ec_step_t*)_data;
		int i;
		for (i = 0; i < data->n_seqs; ++i) {
			bseq1_t *s = &data->seqs[i];
			int is_fq = (s->qual && !es->opt->no_qual);
			if (!es->opt->filter_mode) {
				if (es->opt->discard && (s->aux&7)) goto bfc_ec_cb_free;
				printf("%c%s", is_fq? '@' : '>', s->name);
				if (!s->comment) {
					printf("\tec:Z:%d", s->aux&7);
					if ((s->aux&7) == 0)
						printf("_%d:%d_%d_%d:%d_%d", s->aux2>>10, s->aux2&0xff, s->aux>>3&1, s->aux>>18&0x3fff, s->aux>>4&0x3fff, s->aux2>>8&3);
				} else printf("\t%s", s->comment);
			} else {
				if (s->aux) goto bfc_ec_cb_free;
				printf("%c%s", is_fq? '@' : '>', s->name);
				if (s->comment) printf("\t%s", s->comment);
			}
			putchar('\n'); puts(s->seq);
			if (is_fq) { puts("+"); puts(s->qual); }
bfc_ec_cb_free:
			free(s->seq); free(s->qual); free(s->comment); free(s->name);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

void bfc_correct(const char *fn, const bfc_opt_t *opt, const void *ptr)
{
	ec_shared_t es;
	memset(&es, 0, sizeof(ec_shared_t));
	es.opt = opt;
	if (bfc_verbose >= 3)
		fprintf(stderr, "[M::%s @%.1f*%.1f%%] Starting...\n", __func__, realtime()-bfc_real_time,
				100.*cputime()/(realtime()-bfc_real_time+1e-6));
	if (!opt->filter_mode) {
		int i, mode;
		const bfc_ch_t *ch = (const bfc_ch_t*)ptr;
		uint64_t hist[256], hist_high[64];

		mode = bfc_ch_hist(ch, hist, hist_high);
		if (bfc_verbose >= 4) {
			for (i = 0; i < 256; ++i)
				if (i < 64) fprintf(stderr, "[M::%s] %3d : %llu : %llu\n", __func__, i, (long long)hist[i], (long long)hist_high[i]);
				else fprintf(stderr, "[M::%s] %3d : %llu\n", __func__, i, (long long)hist[i]);
		}

		es.e = calloc(opt->n_threads, sizeof(void*));
		for (i = 0; i < opt->n_threads; ++i)
			es.e[i] = ec1buf_init(opt, ch), es.e[i]->mode = mode;
		es.ks = bseq_open(fn);
		kt_pipeline(opt->no_mt_io? 1 : 2, bfc_ec_cb, &es, 3);
		bseq_close(es.ks);
		for (i = 0; i < opt->n_threads; ++i)
			ec1buf_destroy(es.e[i]);
		free(es.e);
	} else {
		es.bf = (const bfc_bf_t*)ptr;
		es.ks = bseq_open(fn);
		kt_pipeline(opt->no_mt_io? 1 : 2, bfc_ec_cb, &es, 3);
		bseq_close(es.ks);
	}
}
