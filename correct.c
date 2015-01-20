#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include "bfc.h"

static inline int bfc_trusted(const bfc_bf_t *bf, int k, const bfc_kmer_t *x)
{
	uint64_t hash, y[2];
	hash = bfc_kmer_hash(k, x->x, y);
	return (bfc_bf_get(bf, hash) == bf->n_hashes);
}

/**************************
 * Sequence struct for ec *
 **************************/

#include "kvec.h"

typedef struct { // NOTE: unaligned memory
	uint8_t b:3, ob:3, q:1, oq:1;
	uint8_t ec:1, absent:1, dummy6:6;
	uint8_t cov:6, solid_end:1, conflict:1;
	uint8_t dummy;
	int i;
} ecbase_t;

typedef kvec_t(ecbase_t) ecseq_t;

static int bfc_seq_conv(const char *s, const char *q, int qthres, ecseq_t *seq)
{
	int i, l;
	l = strlen(s);
	kv_resize(ecbase_t, *seq, l);
	seq->n = l;
	for (i = 0; i < l; ++i) {
		ecbase_t *c = &seq->a[i];
		c->b = c->ob = seq_nt6_table[(int)s[i]] - 1;
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

int bfc_ec_greedy_k(int k, const bfc_kmer_t *x, const bfc_bf_t *bf)
{
	int i, j, n_trusted = 0, ret = -1;
	for (i = 0; i < k; ++i) {
		int c = (x->x[1]>>i&1)<<1 | (x->x[0]>>i&1);
		for (j = 0; j < 4; ++j) {
			bfc_kmer_t y = *x;
			if (j == c) continue;
			bfc_kmer_change(k, y.x, i, j);
			if (!bfc_trusted(bf, k, &y)) continue;
			ret = i<<2 | j;
			if (++n_trusted > 1) return -1;
		}
	}
	return ret;
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

void bfc_ec_kcov(int k, int min_occ, ecseq_t *s, const bfc_bf_t *bf)
{
	int i, l, j;
	bfc_kmer_t x = bfc_kmer_null;
	for (i = l = 0; i < s->n; ++i) {
		ecbase_t *c = &s->a[i];
		c->solid_end = c->cov = 0;
		if (c->b < 4) {
			bfc_kmer_append(k, x.x, c->b);
			if (++l >= k) {
				if (bfc_trusted(bf, k, &x)) {
					c->solid_end = 1;
					for (j = i - k + 1; j <= i; ++j) ++s->a[j].cov;
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

typedef struct {
	uint8_t ec:1, ec_high:1, absent:1, dummy:1, b:4;
} bfc_penalty_t;

typedef struct {
	int n_ec, n_ec_high, n_absent;
} ecstats_t;

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
} ecstack1_t;

typedef struct {
	const bfc_opt_t *opt;
	const bfc_bf_t *bf;
	kvec_t(echeap1_t) heap;
	kvec_t(ecstack1_t) stack;
	ecseq_t seq, tmp, ec[2];
} bfc_ec1buf_t;

#define heap_lt(a, b) ((a).tot_pen > (b).tot_pen)
KSORT_INIT(ec, echeap1_t, heap_lt)

static bfc_ec1buf_t *ec1buf_init(const bfc_opt_t *opt, const bfc_bf_t *bf)
{
	bfc_ec1buf_t *e;
	e = calloc(1, sizeof(bfc_ec1buf_t));
	e->opt = opt, e->bf = bf;
	return e;
}

static void ec1buf_destroy(bfc_ec1buf_t *e)
{	
	free(e->heap.a); free(e->stack.a); free(e->seq.a); free(e->tmp.a); free(e->ec[0].a); free(e->ec[1].a);
	free(e);
}

#define weighted_penalty(o, p) ((o)->w_ec * (p).ec + (o)->w_ec_high * (p).ec_high + (o)->w_absent * (p).absent)

static void buf_update(bfc_ec1buf_t *e, const echeap1_t *prev, bfc_penalty_t pen)
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

static void buf_backtrack(ecstack1_t *s, int end, const ecseq_t *seq, ecseq_t *path)
{
	int i;
	kv_resize(ecbase_t, *path, seq->n);
	path->n = seq->n;
	while (end >= 0) {
		i = s[end].i;
		path->a[i].b = s[end].b;
		path->a[i].ec = s[end].pen.ec;
		path->a[i].absent = s[end].pen.absent;
		end = s[end].parent;
	}
}

static int bfc_ec1dir(bfc_ec1buf_t *e, const ecseq_t *seq, ecseq_t *ec, int start, int end)
{
	echeap1_t z;
	int i, l, path[BFC_MAX_PATHS], n_paths = 0, n_failures = 0, min_path = -1, min_path_pen = INT_MAX;
	assert(end <= seq->n && end - start >= e->opt->k);
	if (bfc_verbose >= 4) fprintf(stderr, "* bfc_ec1dir(): len:%ld start:%d end:%d\n", seq->n, start, end);
	e->heap.n = e->stack.n = 0;
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
		if (e->heap.n == 0) { // may happen when there is an uncorrectable "N"
			if (n_paths) break;
			else return -2;
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
			int b, os = -1, fixed = 0, other_ext = 0, n_added = 0;
			bfc_penalty_t added[4];
			// test if the read extension alone is enough
			if (z.i > end) fixed = 1;
			if (c && c->b < 4) { // A, C, G or T
				bfc_kmer_t x = z.x;
				bfc_kmer_append(e->opt->k, x.x, c->b);
				os = bfc_trusted(e->bf, e->opt->k, &x);
				if (c->q && os && c->cov >= e->opt->min_cov + 1) fixed = 1;
				if (bfc_verbose >= 4)
					fprintf(stderr, "     Original base:%c qual:%d fixed:%d trusted:%d", "ACGTN"[c->b], c->q, fixed, os);
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
					s = bfc_trusted(e->bf, e->opt->k, &x);
					if (bfc_verbose >= 4 && s)
						fprintf(stderr, "     Alternative k-mer: (%c,%d)\n", "ACGTN"[b], s);
					if (s == 0) continue; // not solid
					pen.ec = c && c->b < 4? 1 : 0;
					pen.ec_high = pen.ec? c->oq : 0;
					pen.absent = 0;
					pen.b = b;
					added[n_added++] = pen;
					++other_ext;
				} else {
					pen.ec = pen.ec_high = 0;
					pen.absent = (os == 0);
					pen.b = b;
					added[n_added++] = pen;
				}
			} // ~for(b)
			if (fixed == 0 && other_ext == 0) ++n_failures;
			if (n_failures > seq->n * 2) {
				if (bfc_verbose >= 4) fprintf(stderr, "  !! too many unsuccessful attempts\n");
				break;
			}
			if (c || n_added == 1) {
				if (n_added > 1 && e->heap.n > e->opt->max_heap) { // to prevent heap explosion
					int min_b = -1, min = INT_MAX;
					for (b = 0; b < n_added; ++b) {
						int t = weighted_penalty(e->opt, added[b]);
						if (min > t) min = t, min_b = b;
					}
					buf_update(e, &z, added[min_b]);
				} else {
					for (b = 0; b < n_added; ++b)
						buf_update(e, &z, added[b]);
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
	if (n_paths == 0) return -3;
	assert(min_path >= 0 && min_path < n_paths && e->stack.a[path[min_path]].tot_pen == min_path_pen);
	buf_backtrack(e->stack.a, path[min_path], seq, ec);
	for (i = 0; i < ec->n; ++i) // mask out uncorrected regions
		if (i < start + e->opt->k || i >= end) ec->a[i].b = 4;
	if (bfc_verbose >= 4) {
		fprintf(stderr, "* %d path(s); lowest penalty: %d\n  ", n_paths, min_path_pen);
		for (i = 0; i < ec->n; ++i) fputc((seq->a[i].b == ec->a[i].b? "ACGTN":"acgtn")[ec->a[i].b], stderr);
		fputc('\n', stderr);
	}
	return 0;
}

typedef struct {
	uint32_t failed:1, n_ec:16, n_ec_high:15;
} ecstat_t;

ecstat_t bfc_ec1(bfc_ec1buf_t *e, char *seq, char *qual)
{
	int i, start = 0, end = 0, n_n = 0;
	uint64_t r;
	ecstat_t s;

	s.failed = 1, s.n_ec = 0, s.n_ec_high = 0;
	bfc_seq_conv(seq, qual, e->opt->q, &e->seq);
	for (i = 0; i < e->seq.n; ++i)
		if (e->seq.a[i].ob > 3) ++n_n;
	if (n_n > e->seq.n * .05) return s;
	bfc_ec_kcov(e->opt->k, e->opt->min_cov, &e->seq, e->bf);
	r = bfc_ec_best_island(e->opt->k, &e->seq);
	if (r == 0) { // no solid k-mer
		bfc_kmer_t x;
		int ec = -1;
		while ((end = bfc_ec_first_kmer(e->opt->k, &e->seq, start, &x)) < e->seq.n) {
			ec = bfc_ec_greedy_k(e->opt->k, &x, e->bf);
			if (ec >= 0) break;
			if (end + (e->opt->k>>1) >= e->seq.n) break;
			start = end - (e->opt->k>>1);
		}
		if (ec >= 0) {
			e->seq.a[end - (ec>>2)].b = ec&3;
			++end; start = end - e->opt->k;
		} else return s; // cannot find a solid k-mer anyway
	} else start = r>>32, end = (uint32_t)r;
	if (bfc_verbose >= 4)
		fprintf(stderr, "* Longest solid island: [%d,%d)\n", start, end);
	if (bfc_ec1dir(e, &e->seq, &e->ec[0], start, e->seq.n) < 0) return s;
	bfc_seq_revcomp(&e->seq);
	if (bfc_ec1dir(e, &e->seq, &e->ec[1], e->seq.n - end, e->seq.n) < 0) return s;
	s.failed = 0;
	bfc_seq_revcomp(&e->ec[1]);
	bfc_seq_revcomp(&e->seq);
	for (i = 0; i < e->seq.n; ++i) {
		ecbase_t *c = &e->seq.a[i];
		if (e->ec[0].a[i].b == e->ec[1].a[i].b)
			c->b = e->ec[0].a[i].b > 3? e->seq.a[i].b : e->ec[0].a[i].b;
		else if (e->ec[1].a[i].b > 3) c->b = e->ec[0].a[i].b;
		else if (e->ec[0].a[i].b > 3) c->b = e->ec[1].a[i].b;
		else c->b = e->seq.a[i].ob, c->conflict = 1;
	}
	for (i = 0; i < e->seq.n; ++i) {
		int is_diff = !(e->seq.a[i].b == e->seq.a[i].ob);
		if (is_diff) {
			++s.n_ec;
			if (e->seq.a[i].q) ++s.n_ec_high;
		}
		seq[i] = (is_diff? "acgtn" : "ACGTN")[e->seq.a[i].b];
		if (qual) qual[i] = is_diff? '+' : "+?"[e->seq.a[i].q];
	}
	if (bfc_verbose >= 4) {
		bfc_ec_kcov(e->opt->k, e->opt->min_cov, &e->seq, e->bf);
		fprintf(stderr, "* failed:%d n_ec:%d n_ec_high:%d\n  ", s.failed, s.n_ec, s.n_ec_high);
		for (i = 0; i < e->seq.n; ++i)
			fputc((e->seq.a[i].b == e->seq.a[i].ob? "ACGTN" : "acgtn")[e->seq.a[i].b], stderr);
		fprintf(stderr, "\n  ");
		for (i = 0; i < e->seq.n; ++i)
			fputc('0' + (int)(10. * e->seq.a[i].cov / e->opt->k + .499), stderr);
		fputc('\n', stderr);
	}
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
	bfc_ec1buf_t **e;
	int64_t n_processed;
} ec_shared_t;

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	ec_shared_t *es;
} ec_step_t;

static void worker_ec(void *_data, long k, int tid)
{
	ec_step_t *data = (ec_step_t*)_data;
	ec_shared_t *es = data->es;
	bseq1_t *s = &data->seqs[k];
	if (bfc_verbose >= 4) fprintf(stderr, "* Processing read '%s'...\n", s->name);
	if (es->opt->kmer_trim) {
		uint64_t max;
		max = max_streak(es->opt->k, es->e[tid]->bf, s);
		if (max>>32 && (double)((max>>32) + es->opt->k) / s->l_seq > es->opt->trim_thres) {
			int start = (uint32_t)max, end = start + (max>>32);
			start -= es->opt->k - 1;
			assert(start >= 0 && end <= s->l_seq);
			if (start > 0) memmove(s->seq, s->seq + start, end - start);
			s->l_seq = end - start;
			s->seq[s->l_seq] = 0;
			if (s->qual) {
				if (start > 0) memmove(s->qual, s->qual + start, s->l_seq);
				s->qual[s->l_seq] = 0;
			}
			s->aux = 0;
		} else s->aux = 1;
	} else {
		ecstat_t st;
		st = bfc_ec1(es->e[tid], s->seq, s->qual);
		s->aux = st.n_ec<<16 | st.n_ec_high<<1 | st.failed;
	}
}

void *bfc_ec_cb(void *shared, int step, void *_data)
{
	ec_shared_t *es = (ec_shared_t*)shared;
	if (step == 0) {
		ec_step_t *ret;
		ret = calloc(1, sizeof(ec_step_t));
		ret->seqs = bseq_read(es->ks, es->opt->chunk_size, &ret->n_seqs);
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
			if (!es->opt->discard || !(s->aux&1)) {
				printf("%c%s\tec:Z:%c_%d_%d\n%s\n", is_fq? '@' : '>', s->name,
					   "TF"[s->aux&1], s->aux>>16&0xffff, s->aux>>1&0x7fff, s->seq);
				if (is_fq) printf("+\n%s\n", s->qual);
			}
			free(s->seq); free(s->qual); free(s->name);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

void bfc_correct(const char *fn, const bfc_opt_t *opt, const bfc_bf_t *bf)
{
	ec_shared_t es;
	int i;

	memset(&es, 0, sizeof(ec_shared_t));
	es.opt = opt;
	if (bfc_verbose >= 3)
		fprintf(stderr, "[M::%s @%.1f*%.1f%%] Starting...\n", __func__, realtime()-bfc_real_time,
				100.*cputime()/(realtime()-bfc_real_time+1e-6));
	es.e = calloc(opt->n_threads, sizeof(void*));
	for (i = 0; i < opt->n_threads; ++i)
		es.e[i] = ec1buf_init(opt, bf);
	es.ks = bseq_open(fn);
	kt_pipeline(opt->no_mt_io? 1 : 2, bfc_ec_cb, &es, 3);
	bseq_close(es.ks);
	for (i = 0; i < opt->n_threads; ++i)
		ec1buf_destroy(es.e[i]);
	free(es.e);
}
