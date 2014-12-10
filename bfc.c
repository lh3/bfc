#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <assert.h>

/******************
 * Hash functions *
 ******************/

// Thomas Wang's integer hash functions. See <https://gist.github.com/lh3/59882d6b96166dfc3d8d> for a snapshot.
static inline uint64_t bfc_hash_64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

/****************
 * Sequence I/O *
 ****************/

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

unsigned char seq_nt6_table[256] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

typedef struct {
	int l_seq;
	char *name, *seq, *qual;
} bseq1_t;

bseq1_t *bseq_read(kseq_t *ks, int chunk_size, int *n_)
{
	int size = 0, m, n;
	bseq1_t *seqs;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		bseq1_t *s;
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
		s->name = strdup(ks->name.s);
		s->seq = strdup(ks->seq.s);
		s->qual = ks->qual.l? strdup(ks->qual.s) : 0;
		s->l_seq = ks->seq.l;
		size += seqs[n++].l_seq;
		if (size >= chunk_size) break;
	}
	*n_ = n;
	return seqs;
}

/**********
 * Timers *
 **********/

#include <sys/resource.h>
#include <sys/time.h>

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

/*****************
 * Configuration *
 *****************/

#include <math.h>

int bfc_verbose = 3;

typedef struct {
	int chunk_size;
	int n_threads;
	int k, q;
	int n_shift, n_hashes;
	int min_cov; // a k-mer is considered solid if the count is no less than this
	int win_multi_ec; // no 2 high-qual corrections or 4 corrections in a window of this size
} bfc_opt_t;

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
	opt->win_multi_ec = 16;
}

void bfc_opt_by_size(bfc_opt_t *opt, long size)
{
	opt->k = (int)(log(size) / log(2) * 1.2);
	if (opt->k > 37) opt->k = 37;
	opt->n_shift = opt->k;
}

/**********************
 * Basic k-mer update *
 **********************/

typedef struct {
	uint64_t x[4];
} bfc_kmer_t;

static bfc_kmer_t bfc_kmer_null = {{0,0,0,0}};

static inline void bfc_kmer_append(int k, uint64_t x[4], int c)
{ // IMPORTANT: 0 <= c < 4
	uint64_t mask = (1ULL<<k) - 1;
	x[0] = (x[0]<<1 | (c&1))  & mask;
	x[1] = (x[1]<<1 | (c>>1)) & mask;
	x[2] = x[2]>>1 | (1ULL^(c&1))<<(k-1);
	x[3] = x[3]>>1 | (1ULL^c>>1) <<(k-1);
}

/* A note on multi-threading

   The bloom filter is always the same regardless of how many threads in use.
   However, the k-mer inserted to the hash table may be slightly different.
   Suppose k-mers A and B are both singletons and that if A is inserted first,
   B is a false positive and gets inserted to the hash table. In the
   multi-threading mode, nonetheless, B may be inserted before A. In this case,
   B is not a false positive any more. This is not a bug. The k-mers put into
   the hash table depends on the order of input.
*/

/************************
 * Blocked Bloom Filter *
 ************************/

#define BFC_BLK_SHIFT  9 // 64 bytes, the size of a cache line
#define BFC_BLK_MASK   ((1<<(BFC_BLK_SHIFT)) - 1)

typedef struct {
	int n_shift, n_hashes;
	uint8_t *b;
} bfc_bf_t;

bfc_bf_t *bfc_bf_init(int n_shift, int n_hashes)
{
	bfc_bf_t *b;
	if (n_shift + BFC_BLK_SHIFT > 64 || n_shift < BFC_BLK_SHIFT) return 0;
	b = calloc(1, sizeof(bfc_bf_t));
	b->n_shift = n_shift;
	b->n_hashes = n_hashes;
	posix_memalign((void**)&b->b, 1<<(BFC_BLK_SHIFT-3), 1ULL<<(n_shift-3));
	bzero(b->b, 1ULL<<(n_shift-3));
	return b;
}

void bfc_bf_destroy(bfc_bf_t *b)
{
	if (b == 0) return;
	free(b->b); free(b);
}

int bfc_bf_insert(bfc_bf_t *b, uint64_t hash)
{
	int x = b->n_shift - BFC_BLK_SHIFT;
	uint64_t y = hash & ((1ULL<<x) - 1);
	int h1 = hash >> x & BFC_BLK_MASK;
	int h2 = hash >> b->n_shift & BFC_BLK_MASK;
	uint8_t *p = &b->b[y<<(BFC_BLK_SHIFT-3)];
	int i, z = h1, cnt = 0;
	if ((h2&31) == 0) h2 = (h2 + 1) & BFC_BLK_MASK; // otherwise we may repeatedly use a few bits
	while (__sync_lock_test_and_set(p, 1)); // lock
	for (i = 0; i < b->n_hashes; z = (z + h2) & BFC_BLK_MASK) {
		uint8_t *q = &p[z>>3], u;
		if (p == q) continue; // don't use the first byte. It is a spin lock.
		u = 1ULL<<(z&7);
		cnt += !!(*q & u);
		*q |= u;
		++i;
	}
	__sync_lock_release(p); // unlock
	return cnt;
}

/**************
 * Hash table *
 **************/

#include "khash.h"

#define _cnt_eq(a, b) ((a)>>14 == (b)>>14)
#define _cnt_hash(a) ((a)>>14)
KHASH_INIT(cnt, uint64_t, char, 0, _cnt_hash, _cnt_eq)
typedef khash_t(cnt) cnthash_t;

#define BFC_CH_KEYBITS 50

#define bfc_flip_cnt(r) (((r) & 0xff) | ((r)>>8&7)<<11 | ((r)>>11&7)<<8)

typedef struct {
	int k;
	cnthash_t **h;
	// private
	int l_pre;
} bfc_ch_t;

bfc_ch_t *bfc_ch_init(int k)
{
	bfc_ch_t *ch;
	int i;
	if (k < 2) return 0;
	ch = calloc(1, sizeof(bfc_ch_t));
	ch->k = k;
	ch->l_pre = k*2 - BFC_CH_KEYBITS; // TODO: this should be improved!!!
	if (ch->l_pre < 1) ch->l_pre = 1;
	if (ch->l_pre > k - 1) ch->l_pre = k - 1;
	ch->h = calloc(1<<ch->l_pre, sizeof(void*));
	for (i = 0; i < 1<<ch->l_pre; ++i)
		ch->h[i] = kh_init(cnt);
	return ch;
}

void bfc_ch_destroy(bfc_ch_t *ch)
{
	int i;
	if (ch == 0) return;
	for (i = 0; i < 1<<ch->l_pre; ++i)
		kh_destroy(cnt, ch->h[i]);
	free(ch->h); free(ch);
}

void bfc_ch_insert(bfc_ch_t *ch, uint64_t x[2], int low_flag)
{
	int absent;
	cnthash_t *h = ch->h[x[0] & ((1ULL<<ch->l_pre) - 1)];
	uint64_t key = (x[0] >> ch->l_pre | x[1] << (ch->k - ch->l_pre)) << 14 | 1;
	khint_t k;
	if (__sync_lock_test_and_set(&h->lock, 1))
		while (__sync_lock_test_and_set(&h->lock, 1))
			while (h->lock); // lock
	k = kh_put(cnt, h, key, &absent);
	if (absent) {
		if (low_flag & 1) kh_key(h, k) |= 1<<8;
		if (low_flag & 2) kh_key(h, k) |= 1<<11;
	} else {
		if ((kh_key(h, k) & 0xff) != 0xff) ++kh_key(h, k);
		if ((low_flag & 1) && (kh_key(h, k) >> 8 & 7) != 7) kh_key(h, k) += 1<<8;
		if ((low_flag & 2) && (kh_key(h, k) >> 11& 7) != 7) kh_key(h, k) += 1<<11;
	}
	__sync_lock_release(&h->lock); // unlock
}

int bfc_ch_get(const bfc_ch_t *ch, const uint64_t x[2])
{
	uint64_t key;
	cnthash_t *h;
	khint_t itr;
	h = ch->h[x[0] & ((1ULL<<ch->l_pre) - 1)];
	key = (x[0] >> ch->l_pre | x[1] << (ch->k - ch->l_pre)) << 14 | 1;
	itr = kh_get(cnt, h, key);
	return itr == kh_end(h)? -1 : kh_key(h, itr) & 0x3fff;
}

uint64_t bfc_ch_count(const bfc_ch_t *ch)
{
	int i;
	uint64_t cnt = 0;
	for (i = 0; i < 1<<ch->l_pre; ++i)
		cnt += kh_size(ch->h[i]);
	return cnt;
}

int bfc_ch_hist(const bfc_ch_t *ch, uint64_t cnt[256])
{
	int i, max_i;
	uint64_t max;
	memset(cnt, 0, 256 * 8);
	for (i = 0; i < 1<<ch->l_pre; ++i) {
		khint_t k;
		cnthash_t *h = ch->h[i];
		for (k = 0; k != kh_end(h); ++k)
			if (kh_exist(h, k))
				++cnt[kh_key(h, k) & 0xff];
	}
	for (i = 0, max = 0; i < 256; ++i) {
		if (cnt[i] > max)
			max = cnt[i], max_i = i;
		if (bfc_verbose >= 3 && i > 0)
			fprintf(stderr, "[M::%s] %3d : %lld\n", __func__, (int)i, (long long)cnt[i]);
	}
	return max_i;
}

int bfc_ch_dump(const bfc_ch_t *ch, const char *fn)
{
	FILE *fp;
	uint32_t t[2];
	int i;
	if ((fp = strcmp(fn, "-")? fopen(fn, "wb") : stdout) == 0) return -1;
	t[0] = ch->k, t[1] = ch->l_pre;
	fwrite(t, 4, 2, fp);
	for (i = 0; i < 1<<ch->l_pre; ++i) {
		cnthash_t *h = ch->h[i];
		khint_t k;
		t[0] = kh_n_buckets(h), t[1] = kh_size(h);
		fwrite(t, 4, 2, fp);
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k))
				fwrite(&kh_key(h, k), 8, 1, fp);
	}
	fprintf(stderr, "[M::%s] dumpped the hash table to file '%s'.\n", __func__, fn);
	fclose(fp);
	return 0;
}

bfc_ch_t *bfc_ch_restore(const char *fn)
{
	FILE *fp;
	uint32_t t[2];
	int i, j, absent;
	bfc_ch_t *ch;

	if ((fp = fopen(fn, "rb")) == 0) return 0;
	fread(t, 4, 2, fp);
	ch = bfc_ch_init(t[0]);
	assert((int)t[1] == ch->l_pre);
	for (i = 0; i < 1<<ch->l_pre; ++i) {
		cnthash_t *h = ch->h[i];
		fread(t, 4, 2, fp);
		kh_resize(cnt, h, t[0]);
		for (j = 0; j < t[1]; ++j) {
			uint64_t key;
			fread(&key, 8, 1, fp);
			kh_put(cnt, h, key, &absent);
			assert(absent);
		}
	}
	fclose(fp);
	fprintf(stderr, "[M::%s] restored the hash table to file '%s'.\n", __func__, fn);
	return ch;
}

/**************
 * Cache hash *
 **************/

typedef struct {
	uint64_t h0:56, ch:8;
	uint64_t h1:56, cl:6, absent:2;
} bfc_kcelem_t;

#define kc_hash(a) ((a).h0 + (a).h1)
#define kc_equal(a, b) ((a).h0 == (b).h0 && (a).h1 == (b).h1)
KHASH_INIT(kc, bfc_kcelem_t, char, 0, kc_hash, kc_equal)
typedef khash_t(kc) kchash_t;

int bfc_kc_get(const bfc_ch_t *ch, kchash_t *kc, const bfc_kmer_t *z)
{
	int r, flipped = !(z->x[0] + z->x[1] < z->x[2] + z->x[3]);
	uint64_t x[2], mask = (1ULL<<ch->k) - 1;
	x[0] = bfc_hash_64(z->x[flipped<<1|0], mask);
	x[1] = bfc_hash_64(z->x[flipped<<1|1], mask);
	if (kc) {
		khint_t k;
		int absent;
		bfc_kcelem_t key, *p;
		key.h0 = x[0], key.h1 = x[1];
		k = kh_put(kc, kc, key, &absent);
		p = &kh_key(kc, k);
		if (absent) {
			r = bfc_ch_get(ch, x);
			if (r >= 0) p->ch = r&0xff, p->cl = r>>8&0x3f, p->absent = 0;
			else p->ch = p->cl = 0, p->absent = 1;
		} else r = p->absent? -1 : p->cl<<8 | p->ch;
	} else r = bfc_ch_get(ch, x);
	if (r >= 0 && flipped)
		r = (r & 0xff) | (r>>8&7)<<11 | (r>>11&7)<<8;
	return r;
}

void bfc_kc_print_kcov(const bfc_ch_t *ch, kchash_t *kc, const char *seq)
{
	int len, i, l, r, c;
	len = strlen(seq);
	bfc_kmer_t x = bfc_kmer_null;
	for (i = l = 0; i < len; ++i) {
		if ((c = seq_nt6_table[(uint8_t)seq[i]] - 1) < 4) {
			bfc_kmer_append(ch->k, x.x, c);
			if (++l >= ch->k && (r = bfc_kc_get(ch, kc, &x)) >= 0)
				fprintf(stderr, "%d\t%d\t%d\t%d\t%d\n", i - ch->k + 1, i + 1, r&0xff, r>>11&7, r>>8&7);
		} else l = 0, x = bfc_kmer_null;
	}
}

/*******************************
 * Other routines for counting *
 *******************************/

static double bfc_real_time;

typedef struct {
	const bfc_opt_t *opt;
	uint64_t mask;
	kseq_t *ks;
	bfc_bf_t *bf;
	bfc_ch_t *ch;
} bfc_cntaux_t;

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	bfc_cntaux_t *aux;
} bfc_count_data_t;

void bfc_kmer_insert(bfc_cntaux_t *aux, const bfc_kmer_t *x, int low_flag)
{
	int k = aux->opt->k, ret;
	int t = !(x->x[0] + x->x[1] < x->x[2] + x->x[3]);
	uint64_t mask = (1ULL<<k) - 1, y[2], hash;
	y[0] = bfc_hash_64(x->x[t<<1|0], mask);
	y[1] = bfc_hash_64(x->x[t<<1|1], mask);
	hash = (y[0] ^ y[1]) << k | ((y[0] + y[1]) & mask);
	ret = bfc_bf_insert(aux->bf, hash);
	if (t) low_flag = (low_flag&1)<<1 | (low_flag&2)>>1;
	if (ret == aux->opt->n_hashes)
		bfc_ch_insert(aux->ch, y, low_flag);
}

static void worker_count(void *_data, long k, int tid)
{
	bfc_count_data_t *data = (bfc_count_data_t*)_data;
	bfc_cntaux_t *aux = data->aux;
	bseq1_t *s = &data->seqs[k];
	const bfc_opt_t *o = aux->opt;
	int i, l;
	bfc_kmer_t x = bfc_kmer_null;
	for (i = l = 0; i < s->l_seq; ++i) {
		int c = seq_nt6_table[(uint8_t)s->seq[i]] - 1;
		if (c < 4) {
			bfc_kmer_append(o->k, x.x, c);
			if (++l >= o->k) {
				int low_flag = 0;
				if (s->qual) {
					if (s->qual[i] - 33 < o->q) low_flag |= 1;
					if (s->qual[i - o->k + 1] - 33 < o->q) low_flag |= 2;
				}
				bfc_kmer_insert(aux, &x, low_flag);
			}
		} else l = 0, x = bfc_kmer_null;
	}
}

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

void *bfc_count_cb(void *shared, int step, void *_data)
{
	bfc_cntaux_t *aux = (bfc_cntaux_t*)shared;
	if (step == 0) {
		bfc_count_data_t *ret;
		ret = calloc(1, sizeof(bfc_count_data_t));
		ret->seqs = bseq_read(aux->ks, aux->opt->chunk_size, &ret->n_seqs);
		ret->aux = aux;
		fprintf(stderr, "[M::%s] read %d sequences\n", __func__, ret->n_seqs);
		if (ret->seqs) return ret;
		else free(ret);
	} else if (step == 1) {
		int i;
		bfc_count_data_t *data = (bfc_count_data_t*)_data;
		kt_for(aux->opt->n_threads, worker_count, data, data->n_seqs);
		fprintf(stderr, "[M::%s] processed %d sequences (CPU/real time: %.3f/%.3f secs; # distinct k-mers: %ld)\n",
				__func__, data->n_seqs, cputime(), realtime() - bfc_real_time, (long)bfc_ch_count(aux->ch));
		for (i = 0; i < data->n_seqs; ++i) {
			bseq1_t *s = &data->seqs[i];
			free(s->seq); free(s->qual); free(s->name);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

/**************************
 * Sequence struct for ec *
 **************************/

#include "kvec.h"

typedef struct { // NOTE: unaligned memory
	uint8_t b:3, ob:3, q:1, oq:1;
	uint8_t ec:1, absent:1, diff:6;
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

/********************
 * Correct one read *
 ********************/

#include "ksort.h"

#define BFC_MAX_PATHS 8
#define BFC_MAX_PDIFF 3

typedef struct {
	uint8_t ec:1, ec_high:1, absent:1;
} bfc_penalty_t;

typedef struct {
	int n_ec, n_ec_high, n_absent;
} ecstats_t;

typedef struct {
	int tot_pen;
	int i; // base position
	int k; // position in the stack
	int ecpos_high; // position of the last high-qual correction
	uint64_t ecpos4; // positions of last 4 corrections
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
	kchash_t *kc;
	kvec_t(echeap1_t) heap;
	kvec_t(ecstack1_t) stack;
	ecseq_t ori_seq, tmp, ec[2];
} bfc_ec1buf_t;

#define heap_lt(a, b) ((a).tot_pen > (b).tot_pen)
KSORT_INIT(ec, echeap1_t, heap_lt)

static bfc_ec1buf_t *ec1buf_init(const bfc_opt_t *opt, bfc_ch_t *ch)
{
	bfc_ec1buf_t *e;
	e = calloc(1, sizeof(bfc_ec1buf_t));
	e->opt = opt, e->ch = ch;
	e->kc = kh_init(kc);
	return e;
}

static void ec1buf_destroy(bfc_ec1buf_t *e)
{	
	kh_destroy(kc, e->kc);
	free(e->heap.a); free(e->stack.a); free(e->tmp.a); free(e->ec[0].a); free(e->ec[1].a);
}

static void buf_update(bfc_ec1buf_t *e, const echeap1_t *prev, int b, bfc_penalty_t pen)
{
	ecstack1_t *q;
	echeap1_t *r;
	// update stack
	kv_pushp(ecstack1_t, e->stack, &q);
	q->parent = prev->k;
	q->i = prev->i;
	q->b = b;
	q->pen = pen;
	q->tot_pen = prev->tot_pen + (int)pen.ec + pen.ec_high + pen.absent;
	// update heap
	kv_pushp(echeap1_t, e->heap, &r);
	r->i = prev->i + 1;
	r->k = e->stack.n - 1;
	r->x = prev->x;
	if (pen.ec_high) r->ecpos_high = prev->i;
	if (pen.ec || pen.ec_high) r->ecpos4 = r->ecpos4<<16 | prev->i;
	r->tot_pen = q->tot_pen;
	bfc_kmer_append(e->opt->k, r->x.x, b);
	if (bfc_verbose >= 4)
		fprintf(stderr, "     <= base:%c penalty:%d\n", pen.ec? "acgtn"[b] : "ACGTN"[b], r->tot_pen);
	ks_heapup_ec(e->heap.n, e->heap.a);
}

static void buf_backtrack(ecstack1_t *s, int end, const ecseq_t *ori_seq, ecseq_t *path)
{
	int i;
	kv_resize(ecbase_t, *path, ori_seq->n);
	path->n = ori_seq->n;
	for (i = 0; i < path->n; ++i) path->a[i] = ori_seq->a[i];
	while (end >= 0) {
		i = s[end].i;
		path->a[i].b = s[end].b;
		path->a[i].ec = s[end].pen.ec;
		path->a[i].absent = s[end].pen.absent;
		end = s[end].parent;
	}
}

static void adjust_min_diff(int diff, ecseq_t *opt, const ecseq_t *sub)
{
	int i;
	diff = diff < 63? diff : 63;
	for (i = 0; i < opt->n; ++i)
		if (opt->a[i].b != sub->a[i].b && opt->a[i].diff < diff)
			opt->a[i].diff = diff;
}

static int bfc_ec1dir(bfc_ec1buf_t *e, const ecseq_t *seq, ecseq_t *ec)
{
	echeap1_t z;
	int i, l, path[BFC_MAX_PATHS], n_paths = 0, n_failures = 0;
	assert(seq->n < 0x10000);
	if (bfc_verbose >= 4) fprintf(stderr, "o bfc_ec1dir(): len:%ld\n", seq->n);
	e->heap.n = e->stack.n = 0;
	memset(&z, 0, sizeof(echeap1_t));
	kv_resize(ecbase_t, *ec, seq->n);
	for (i = 0; i < seq->n; ++i) ec->a[i] = seq->a[i];
	for (z.i = l = 0; z.i < seq->n; ++z.i) {
		int c = seq->a[z.i].ob;
		if (c < 4) {
			if (++l == e->opt->k) break;
			bfc_kmer_append(e->opt->k, z.x.x, c);
		} else l = 0, z.x = bfc_kmer_null;
	}
	if (z.i == seq->n) return -1;
	z.k = -1; z.ecpos_high = -1;
	kv_push(echeap1_t, e->heap, z);
	// exhaustive error correction
	while (1) {
		if (e->heap.n == 0) return -2; // may happen when there is an uncorrectable "N"
		z = e->heap.a[0];
		e->heap.a[0] = kv_pop(e->heap);
		ks_heapdown_ec(0, e->heap.n, e->heap.a);
		if (bfc_verbose >= 4)
			fprintf(stderr, "  => pos:%d stack_size:%ld heap_size:%ld penalty:%d\n",
					z.i, e->stack.n, e->heap.n, z.tot_pen);
		if (n_paths && z.tot_pen > e->stack.a[path[0]].tot_pen + BFC_MAX_PDIFF) break;
		if (z.i >= seq->n) { // reach to the end
			if (bfc_verbose >= 4) fprintf(stderr, "  ** reached the end\n");
			path[n_paths++] = z.k;
			if (n_paths == BFC_MAX_PATHS) break;
		} else {
			ecbase_t *c = &seq->a[z.i];
			int b, os = -1, no_others = 0, other_ext = 0;
			// test if the read extension alone is enough
			if (c->ob < 4) { // A, C, G or T
				bfc_kmer_t x = z.x;
				bfc_kmer_append(e->opt->k, x.x, c->ob);
				os = bfc_kc_get(e->ch, e->kc, &x);
				if (os >= 0 && (os&0xff) >= e->opt->min_cov && c->oq) no_others = 1;
				if (bfc_verbose >= 4) {
					fprintf(stderr, "     Original k-mer count: %c,", "ACGTN"[c->ob]);
					if (os >= 0) fprintf(stderr, "%d:%d:%d\n", os&0xff, os>>11&7, os>>8&7);
					else fprintf(stderr, "-1:-1:-1\n");
				}
			}
			// extension
			for (b = 0; b < 4; ++b) {
				bfc_penalty_t pen;
				if (no_others && b != c->ob) continue;
				if (b != c->ob) {
					int s;
					bfc_kmer_t x = z.x;
					if (z.ecpos_high >= 0 && z.i - z.ecpos_high < e->opt->win_multi_ec) continue;
					if (z.i - (z.ecpos4>>48) < e->opt->win_multi_ec) continue;
					bfc_kmer_append(e->opt->k, x.x, b);
					s = bfc_kc_get(e->ch, e->kc, &x);
					if (bfc_verbose >= 4 && s >= 0)
						fprintf(stderr, "     Alternative k-mer count: %c,%d:%d:%d\n", "ACGTN"[b], s&0xff, s>>11&7, s>>8&7);
					if (s < 0 || (s&0xff) < e->opt->min_cov) continue; // not solid
					if (os >= 0 && (s&0xff) - (os&0xff) < 2) continue; // not sufficiently good
					pen.ec = 1, pen.ec_high = c->oq;
					pen.absent = 0;
					buf_update(e, &z, b, pen);
					++other_ext;
				} else {
					pen.ec = pen.ec_high = 0;
					pen.absent = (os < 0 || (os&0xff) < e->opt->min_cov);
					buf_update(e, &z, b, pen);
				}
			} // ~for(b)
			if (no_others == 0 && other_ext == 0) ++n_failures;
			if (n_failures > seq->n) break;
		} // ~else
	} // ~while(1)
	// backtrack
	if (n_paths == 0) return -3;
	buf_backtrack(e->stack.a, path[0], seq, ec);
	if (bfc_verbose >= 4) fprintf(stderr, "o %d path(s)\n", n_paths);
	for (i = 0; i < ec->n; ++i) ec->a[i].diff = 63;
	for (i = 1; i < n_paths; ++i) {
		int diff = e->stack.a[path[i]].tot_pen - e->stack.a[path[0]].tot_pen;
		buf_backtrack(e->stack.a, path[i], seq, &e->tmp);
		adjust_min_diff(diff, ec, &e->tmp);
	}
	return 0;
}

void bfc_ec1(bfc_ec1buf_t *e, char *seq, char *qual)
{
	int i, ret[2];
	bfc_seq_conv(seq, qual, e->opt->q, &e->ori_seq);
	ret[0] = bfc_ec1dir(e, &e->ori_seq, &e->ec[0]);
	bfc_seq_revcomp(&e->ori_seq);
	ret[1] = bfc_ec1dir(e, &e->ori_seq, &e->ec[1]);
	bfc_seq_revcomp(&e->ec[1]);
	bfc_seq_revcomp(&e->ori_seq);
	for (i = 0; i < e->ori_seq.n; ++i) {
		if (e->ec[0].a[i].b == e->ec[1].a[i].b)
			e->ori_seq.a[i].b = e->ec[0].a[i].b;
		else e->ori_seq.a[i].b = e->ori_seq.a[i].ob;
	}
	for (i = 0; i < e->ori_seq.n; ++i)
		seq[i] = (e->ori_seq.a[i].b == e->ori_seq.a[i].ob? "ACGTN" : "acgtn")[e->ori_seq.a[i].b];
}

/********************
 * Error correction *
 ********************/

typedef struct {
	const bfc_opt_t *opt;
	kseq_t *ks;
	bfc_ec1buf_t **e;
	int64_t n_processed;
	int mode;
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
	bfc_ec1(aux->e[tid], s->seq, s->qual);
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
			printf("%c%s\n%s\n", s->qual? '@' : '>', s->name, s->seq);
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

int main(int argc, char *argv[])
{
	gzFile fp;
	bfc_opt_t opt;
	bfc_cntaux_t caux;
	int i, c, mode;
	int no_mt_io = 0, no_ec = 0;
	char *in_hash = 0, *out_hash = 0, *str_kcov = 0;
	uint64_t hist[256];

	bfc_real_time = realtime();
	bfc_opt_init(&opt);
	caux.opt = &opt;
	while ((c = getopt(argc, argv, "v:Ed:k:s:b:L:t:C:h:q:Jr:")) >= 0) {
		if (c == 'k') opt.k = atoi(optarg);
		else if (c == 'C') str_kcov = optarg;
		else if (c == 'd') out_hash = optarg;
		else if (c == 'r') in_hash = optarg;
		else if (c == 'q') opt.q = atoi(optarg);
		else if (c == 'b') opt.n_shift = atoi(optarg);
		else if (c == 't') opt.n_threads = atoi(optarg);
		else if (c == 'h') opt.n_hashes = atoi(optarg);
		else if (c == 'J') no_mt_io = 1; // for debugging kt_pipeline()
		else if (c == 'E') no_ec = 1;
		else if (c == 'v') bfc_verbose = atoi(optarg);
		else if (c == 'L' || c == 's') {
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
	caux.mask = (1ULL<<opt.k) - 1;

	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bfc [options] <in.fq>\n\n");
		fprintf(stderr, "Options: -s FLOAT     approx genome size (k/m/g allowed; change -k and -b) [unset]\n");
		fprintf(stderr, "         -k INT       k-mer length [%d]\n", opt.k);
		fprintf(stderr, "         -t INT       number of threads [%d]\n", opt.n_threads);
		fprintf(stderr, "         -b INT       set Bloom Filter size to pow(2,INT) bits [%d]\n", opt.n_shift);
		fprintf(stderr, "         -h INT       use INT hash functions for Bloom Filter [%d]\n", opt.n_hashes);
		fprintf(stderr, "         -d FILE      dump hash table to FILE [null]\n");
		fprintf(stderr, "         -r FILE      restore hash table from FILE [null]\n");
		fprintf(stderr, "         -E           skip error correction\n");
		fprintf(stderr, "\n");
		return 1;
	}

	if (in_hash == 0) {
		caux.bf = bfc_bf_init(opt.n_shift, opt.n_hashes);
		caux.ch = bfc_ch_init(opt.k);

		fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
		caux.ks = kseq_init(fp);
		kt_pipeline(no_mt_io? 1 : 2, bfc_count_cb, &caux, 2);
		kseq_destroy(caux.ks);
		gzclose(fp);

		bfc_bf_destroy(caux.bf);
		caux.bf = 0;
	} else caux.ch = bfc_ch_restore(in_hash);

	mode = bfc_ch_hist(caux.ch, hist);

	if (str_kcov) bfc_kc_print_kcov(caux.ch, 0, str_kcov);
	if (out_hash) bfc_ch_dump(caux.ch, out_hash);

	if (!no_ec) {
		bfc_ecaux_t eaux;
		eaux.opt = &opt;
		eaux.mode = mode;
		eaux.e = calloc(opt.n_threads, sizeof(void*));
		for (i = 0; i < opt.n_threads; ++i)
			eaux.e[i] = ec1buf_init(&opt, caux.ch);
		fp = gzopen(argv[optind], "r");
		eaux.ks = kseq_init(fp);
		kt_pipeline(no_mt_io? 1 : 2, bfc_ec_cb, &eaux, 3);
		kseq_destroy(eaux.ks);
		gzclose(fp);
		for (i = 0; i < opt.n_threads; ++i)
			ec1buf_destroy(eaux.e[i]);
		free(eaux.e);
	}

	bfc_ch_destroy(caux.ch);

	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - bfc_real_time, cputime());
	return 0;
}
