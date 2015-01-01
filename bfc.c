#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include "bfc.h"

#define BFC_VERSION "r82"

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
double cputime(void);
double realtime(void);

int bfc_verbose = 3;
double bfc_real_time;
bfc_kmer_t bfc_kmer_null = {{0,0,0,0}};

/*****************
 * Configuration *
 *****************/

#include <math.h>

void bfc_opt_init(bfc_opt_t *opt)
{
	memset(opt, 0, sizeof(bfc_opt_t));
	opt->chunk_size = 100000000;
	opt->n_threads = 1;
	opt->k = 33;
	opt->q = 20;
	opt->n_shift = 33;
	opt->n_hashes = 4;

	opt->min_cov = 3;
	opt->win_multi_ec = 10;
	opt->max_end_ext = 5;

	opt->w_ec = 1;
	opt->w_ec_high = 7;
	opt->w_absent = 3;
	opt->w_absent_high = 1;
	opt->max_path_diff = 15;
}

void bfc_opt_by_size(bfc_opt_t *opt, long size)
{
	opt->k = (int)(log(size) / log(2) * 1.2);
	if (opt->k > 37) opt->k = 37;
	opt->n_shift = opt->k;
}

/**************************
 * Sequence struct for ec *
 **************************/

#include "kvec.h"

typedef struct { // NOTE: unaligned memory
	uint8_t b:3, ob:3, q:1, oq:1;
	uint8_t ec:1, absent:1, dummy6:6;
	uint16_t lcov:6, hcov:6, solid_end:1, high_end:1, conflict:1, dummy:1;
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
 * independent ec routines *
 ***************************/

int bfc_ec_greedy_k(int k, int mode, const bfc_kmer_t *x, const bfc_ch_t *ch, bfc_kc_t *kc)
{
	int i, j, max = 0, max_ec = -1, max2 = 0;
	for (i = 0; i < k; ++i) {
		int c = (x->x[1]>>i&1)<<1 | (x->x[0]>>i&1);
		for (j = 0; j < 4; ++j) {
			bfc_kmer_t y = *x;
			int ret;
			if (j == c) continue;
			bfc_kmer_change(k, y.x, i, j);
			ret = bfc_kc_get(ch, kc, &y);
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

void bfc_ec_kcov(int k, int min_occ, ecseq_t *s, const bfc_ch_t *ch, bfc_kc_t *kc)
{
	int i, l, r, j;
	bfc_kmer_t x = bfc_kmer_null;
	for (i = l = 0; i < s->n; ++i) {
		ecbase_t *c = &s->a[i];
		c->high_end = c->solid_end = c->lcov = c->hcov = 0;
		if (c->b < 4) {
			bfc_kmer_append(k, x.x, c->b);
			if (++l >= k) {
				if ((r = bfc_kc_get(ch, kc, &x)) >= 0) {
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

#define BFC_MAX_PATHS 8
#define BFC_EC_HIST 5

typedef struct {
	uint8_t ec:1, ec_high:1, absent:1, absent_high:1, b:4;
} bfc_penalty_t;

typedef struct {
	int n_ec, n_ec_high, n_absent;
} ecstats_t;

typedef struct {
	int tot_pen;
	int i; // base position
	int k; // position in the stack
	int ecpos_high; // position of the last high-qual correction
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
	const bfc_ch_t *ch;
	bfc_kc_t *kc;
	kvec_t(echeap1_t) heap;
	kvec_t(ecstack1_t) stack;
	ecseq_t seq, tmp, ec[2];
	int mode;
} bfc_ec1buf_t;

#define heap_lt(a, b) ((a).tot_pen > (b).tot_pen)
KSORT_INIT(ec, echeap1_t, heap_lt)

static bfc_ec1buf_t *ec1buf_init(const bfc_opt_t *opt, bfc_ch_t *ch)
{
	bfc_ec1buf_t *e;
	e = calloc(1, sizeof(bfc_ec1buf_t));
	e->opt = opt, e->ch = ch;
	e->kc = bfc_kc_init();
	return e;
}

static void ec1buf_destroy(bfc_ec1buf_t *e)
{	
	bfc_kc_destroy(e->kc);
	free(e->heap.a); free(e->stack.a); free(e->seq.a); free(e->tmp.a); free(e->ec[0].a); free(e->ec[1].a);
	free(e);
}

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
	q->tot_pen = prev->tot_pen + o->w_ec * pen.ec + o->w_ec_high * pen.ec_high + o->w_absent * pen.absent + o->w_absent_high * pen.absent_high;
	// update heap
	kv_pushp(echeap1_t, e->heap, &r);
	r->i = prev->i + 1;
	r->k = e->stack.n - 1;
	r->x = prev->x;
	r->ecpos_high = pen.ec_high? prev->i : prev->ecpos_high;
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
	z.k = -1; z.ecpos_high = -1;
	for (i = 0; i < BFC_EC_HIST; ++i) z.ecpos[i] = -1;
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
			fprintf(stderr, "  => pos:%d stack_size:%ld heap_size:%ld penalty:%d last_base:%c ecpos_high:%d ecpos:[%d,%d,%d,%d,%d]\n",
					z.i, e->stack.n, e->heap.n, z.tot_pen, "ACGT"[(z.x.x[1]&1)<<1|(z.x.x[0]&1)], z.ecpos_high,
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
				os = bfc_kc_get(e->ch, e->kc, &x);
				if (c->q && (os&0xff) >= e->opt->min_cov + 1 && c->lcov >= e->opt->min_cov + 1) fixed = 1;
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
						if (c->q && z.ecpos_high >= 0 && z.i - z.ecpos_high < e->opt->win_multi_ec) continue; // no close highQ corrections
						if (z.ecpos[BFC_EC_HIST-1] >= 0 && z.i - z.ecpos[BFC_EC_HIST-1] < e->opt->win_multi_ec) continue; // no clustered corrections
					}
					bfc_kmer_append(e->opt->k, x.x, b);
					s = bfc_kc_get(e->ch, e->kc, &x);
					if (bfc_verbose >= 4 && s >= 0)
						fprintf(stderr, "     Alternative k-mer count: %c,%d:%d\n", "ACGTN"[b], s&0xff, s>>8&0x3f);
					if (s < 0 || (s&0xff) < e->opt->min_cov) continue; // not solid
					//if (os >= 0 && (s&0xff) - (os&0xff) < 2) continue; // not sufficiently better than the read path
					pen.ec = c && c->b < 4? 1 : 0;
					pen.ec_high = pen.ec? c->oq : 0;
					pen.absent = 0;
					pen.absent_high = ((s>>8&0xff) < e->opt->min_cov);
					pen.b = b;
					added[n_added++] = pen;
					++other_ext;
				} else {
					pen.ec = pen.ec_high = 0;
					pen.absent = (os < 0 || (os&0xff) < e->opt->min_cov);
					pen.absent_high = (os < 0 || (os>>8&0xff) < e->opt->min_cov);
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
				for (b = 0; b < n_added; ++b)
					buf_update(e, &z, added[b]);
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
	int i, start = 0, end = 0;
	uint64_t r;
	ecstat_t s;

	bfc_kc_clear(e->kc);
	s.failed = 1, s.n_ec = 0, s.n_ec_high = 0;
	bfc_seq_conv(seq, qual, e->opt->q, &e->seq);
	bfc_ec_kcov(e->opt->k, e->opt->min_cov, &e->seq, e->ch, e->kc);
	r = bfc_ec_best_island(e->opt->k, &e->seq);
	if (r == 0) { // no solid k-mer
		bfc_kmer_t x;
		int ec = -1;
		while ((end = bfc_ec_first_kmer(e->opt->k, &e->seq, start, &x)) < e->seq.n) {
			ec = bfc_ec_greedy_k(e->opt->k, e->mode, &x, e->ch, e->kc);
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
		qual[i] = is_diff? '+' : "+?"[e->seq.a[i].q];
	}
	if (bfc_verbose >= 4) {
		bfc_ec_kcov(e->opt->k, e->opt->min_cov, &e->seq, e->ch, e->kc);
		fprintf(stderr, "* failed:%d n_ec:%d n_ec_high:%d\n  ", s.failed, s.n_ec, s.n_ec_high);
		for (i = 0; i < e->seq.n; ++i)
			fputc((e->seq.a[i].b == e->seq.a[i].ob? "ACGTN" : "acgtn")[e->seq.a[i].b], stderr);
		fprintf(stderr, "\n  ");
		for (i = 0; i < e->seq.n; ++i)
			fputc('0' + (int)(10. * e->seq.a[i].lcov / e->opt->k + .499), stderr);
		fputc('\n', stderr);
	}
	return s;
}

/********************
 * Error correction *
 ********************/

typedef struct {
	const bfc_opt_t *opt;
	bseq_file_t *ks;
	bfc_ec1buf_t **e;
	int64_t n_processed;
} bfc_ecaux_t;

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	bfc_ecaux_t *aux;
} bfc_ec_data_t;

static void worker_ec(void *_data, long k, int tid)
{
	bfc_ec_data_t *data = (bfc_ec_data_t*)_data;
	bfc_ecaux_t *aux = data->aux;
	bseq1_t *s = &data->seqs[k];
	ecstat_t st;
	st = bfc_ec1(aux->e[tid], s->seq, s->qual);
	s->aux = st.n_ec<<16 | st.n_ec_high<<1 | st.failed;
}

void *bfc_ec_cb(void *shared, int step, void *_data)
{
	bfc_ecaux_t *aux = (bfc_ecaux_t*)shared;
	if (step == 0) {
		bfc_ec_data_t *ret;
		ret = calloc(1, sizeof(bfc_ec_data_t));
		ret->seqs = bseq_read(aux->ks, aux->opt->chunk_size, &ret->n_seqs);
		ret->aux = aux;
		fprintf(stderr, "[M::%s] read %d sequences\n", __func__, ret->n_seqs);
		if (ret->seqs) return ret;
		else free(ret);
	} else if (step == 1) {
		bfc_ec_data_t *data = (bfc_ec_data_t*)_data;
		kt_for(aux->opt->n_threads, worker_ec, data, data->n_seqs);
		fprintf(stderr, "[M::%s] processed %d sequences (CPU/real time: %.3f/%.3f secs)\n",
				__func__, data->n_seqs, cputime(), realtime() - bfc_real_time);
		return data;
	} else if (step == 2) {
		bfc_ec_data_t *data = (bfc_ec_data_t*)_data;
		int i;
		for (i = 0; i < data->n_seqs; ++i) {
			bseq1_t *s = &data->seqs[i];
			if (aux->opt->discard && (s->aux&1)) continue;
			printf("%c%s\tec:Z:%c_%d_%d\n%s\n", s->qual? '@' : '>', s->name,
				   "TF"[s->aux&1], s->aux>>16&0xffff, s->aux>>1&0x7fff, s->seq);
			if (s->qual) printf("+\n%s\n", s->qual);
			free(s->seq); free(s->qual); free(s->name);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

/**************
 * Main entry *
 **************/

static void usage(FILE *fp, bfc_opt_t *o)
{
	fprintf(fp, "Usage: bfc [options] <in.fq>\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -s FLOAT     approx genome size (k/m/g allowed; change -k and -b) [unset]\n");
	fprintf(fp, "  -k INT       k-mer length [%d]\n", o->k);
	fprintf(fp, "  -t INT       number of threads [%d]\n", o->n_threads);
	fprintf(fp, "  -b INT       set Bloom filter size to pow(2,INT) bits [%d]\n", o->n_shift);
	fprintf(fp, "  -H INT       use INT hash functions for Bloom filter [%d]\n", o->n_hashes);
	fprintf(fp, "  -d FILE      dump hash table to FILE [null]\n");
	fprintf(fp, "  -E           skip error correction\n");
	fprintf(fp, "  -r FILE      restore hash table from FILE [null]\n");
	fprintf(fp, "  -w INT       no more than %d ec or 2 highQ ec in INT-bp window [%d]\n", BFC_EC_HIST, o->win_multi_ec);
	fprintf(fp, "  -c INT       min k-mer coverage [%d]\n", o->min_cov);
	fprintf(fp, "  -D           discard uncorrectable reads\n");
	fprintf(fp, "  -v           show version number\n");
	fprintf(fp, "  -h           show command line help\n");
}

int main(int argc, char *argv[])
{
	bfc_opt_t opt;
	bfc_bf_t *bf = 0;
	bfc_ch_t *ch = 0;
	int i, c, mode;
	int no_mt_io = 0, no_ec = 0;
	char *in_hash = 0, *out_hash = 0;

	bfc_real_time = realtime();
	bfc_opt_init(&opt);
	while ((c = getopt(argc, argv, "hvV:Ed:k:s:b:L:t:C:H:q:Jr:c:w:D1")) >= 0) {
		if (c == 'k') opt.k = atoi(optarg);
		else if (c == 'd') out_hash = optarg;
		else if (c == 'r') in_hash = optarg;
		else if (c == 'q') opt.q = atoi(optarg);
		else if (c == 'b') opt.n_shift = atoi(optarg);
		else if (c == 't') opt.n_threads = atoi(optarg);
		else if (c == 'H') opt.n_hashes = atoi(optarg);
		else if (c == 'c') opt.min_cov = atoi(optarg);
		else if (c == 'w') opt.win_multi_ec = atoi(optarg);
		else if (c == 'D') opt.discard = 1;
		else if (c == '1') opt.filter_mode = 1;
		else if (c == 'J') no_mt_io = 1; // for debugging kt_pipeline()
		else if (c == 'E') no_ec = 1;
		else if (c == 'V') bfc_verbose = atoi(optarg);
		else if (c == 'h') {
			usage(stdout, &opt);
			return 0;
		} else if (c == 'v') {
			printf("%s\n", BFC_VERSION);
			return 0;
		} else if (c == 'L' || c == 's') {
			double x;
			char *p;
			x = strtod(optarg, &p);
			if (*p == 'G' || *p == 'g') x *= 1e9;
			else if (*p == 'M' || *p == 'm') x *= 1e6;
			else if (*p == 'K' || *p == 'k') x *= 1e3;
			if (c == 's') {
				bfc_opt_by_size(&opt, (long)x + 1);
				fprintf(stderr, "[M::%s] set k to %d\n", __func__, opt.k);
			} else if (c == 'L') opt.chunk_size = (long)x + 1;
		}
	}

	if (optind == argc) {
		usage(stderr, &opt);
		return 1;
	}

	if (opt.filter_mode) bf = (bfc_bf_t*)bfc_count(argv[optind], &opt);
	else if (in_hash) ch = bfc_ch_restore(in_hash);
	else ch = (bfc_ch_t*)bfc_count(argv[optind], &opt);

	if (ch) {
		uint64_t hist[256], hist_high[64];
		mode = bfc_ch_hist(ch, hist, hist_high);
		if (bfc_verbose >= 4) {
			for (i = 0; i < 256; ++i)
				if (i < 64) fprintf(stderr, "[M::%s] %3d : %llu : %llu\n", __func__, i, (long long)hist[i], (long long)hist_high[i]);
				else fprintf(stderr, "[M::%s] %3d : %llu\n", __func__, i, (long long)hist[i]);
		}
		if (out_hash) bfc_ch_dump(ch, out_hash);
		if (!no_ec) {
			bfc_ecaux_t eaux;
			eaux.opt = &opt;
			eaux.e = calloc(opt.n_threads, sizeof(void*));
			for (i = 0; i < opt.n_threads; ++i)
				eaux.e[i] = ec1buf_init(&opt, ch), eaux.e[i]->mode = mode;
			eaux.ks = bseq_open(argv[optind]);
			kt_pipeline(no_mt_io? 1 : 2, bfc_ec_cb, &eaux, 3);
			bseq_close(eaux.ks);
			for (i = 0; i < opt.n_threads; ++i)
				ec1buf_destroy(eaux.e[i]);
			free(eaux.e);
		}
		bfc_ch_destroy(ch);
	} else if (bf) {
		bfc_bf_destroy(bf);
	}

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, BFC_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - bfc_real_time, cputime());
	return 0;
}
