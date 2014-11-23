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

// The inversion of hash_64(). Modified from <https://naml.us/blog/tag/invertible>
static inline uint64_t bfc_hash_64i(uint64_t key, uint64_t mask)
{
	uint64_t tmp;

	// Invert key = key + (key << 31)
	tmp = (key - (key << 31));
	key = (key - (tmp << 31)) & mask;

	// Invert key = key ^ (key >> 28)
	tmp = key ^ key >> 28;
	key = key ^ tmp >> 28;

	// Invert key *= 21
	key = (key * 14933078535860113213ull) & mask;

	// Invert key = key ^ (key >> 14)
	tmp = key ^ key >> 14;
	tmp = key ^ tmp >> 14;
	tmp = key ^ tmp >> 14;
	key = key ^ tmp >> 14;

	// Invert key *= 265
	key = (key * 15244667743933553977ull) & mask;

	// Invert key = key ^ (key >> 24)
	tmp = key ^ key >> 24;
	key = key ^ tmp >> 24;

	// Invert key = (~key) + (key << 21)
	tmp = ~key;
	tmp = ~(key - (tmp << 21));
	tmp = ~(key - (tmp << 21));
	key = ~(key - (tmp << 21)) & mask;

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
	char *seq, *qual;
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

typedef struct {
	int chunk_size;
	int n_threads;
	int k;
	int n_shift, n_hashes;
} bfc_opt_t;

void bfc_opt_init(bfc_opt_t *opt)
{
	memset(opt, 0, sizeof(bfc_opt_t));
	opt->chunk_size = 100000000;
	opt->n_threads = 1;
	opt->k = 33;
	opt->n_shift = 33;
	opt->n_hashes = 4;
}

void bfc_opt_by_size(bfc_opt_t *opt, long size)
{
	opt->k = (int)(log(size) / log(2) * 1.2);
	if (opt->k > 37) opt->k = 37;
	opt->n_shift = opt->k;
}

typedef struct {
	uint64_t x[4];
} bfc_kmer_t;

static bfc_kmer_t bfc_kmer_null = {{0,0,0,0}};

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

uint64_t bfc_ch_lock_cnt = 0;

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

void bfc_ch_insert(bfc_ch_t *ch, uint64_t x[2])
{
	int absent;
	cnthash_t *h = ch->h[x[0] & ((1ULL<<ch->l_pre) - 1)];
	uint64_t key = (x[0] >> ch->l_pre | x[1] << (ch->k - ch->l_pre)) << 14 | 1;
	khint_t k;
	while (__sync_lock_test_and_set(&h->lock, 1)); // lock
	k = kh_put(cnt, h, key, &absent);
	if (!absent && (kh_key(h, k) & 0xff) != 0xff)
		++kh_key(h, k);
	__sync_lock_release(&h->lock); // unlock
}

uint64_t bfc_ch_count(const bfc_ch_t *ch)
{
	int i;
	uint64_t cnt = 0;
	for (i = 0; i < 1<<ch->l_pre; ++i)
		cnt += kh_size(ch->h[i]);
	return cnt;
}

/**********************
 **********************/

typedef struct {
	bfc_opt_t opt;
	uint64_t mask;
	int n_seqs;
	bseq1_t *seqs;
	bfc_bf_t *bf;
	bfc_ch_t *ch;
} bfc_aux_t;

void bfc_kmer_insert(bfc_aux_t *aux, const bfc_kmer_t *x)
{
	int k = aux->opt.k, ret;
	int t = !(x->x[0] + x->x[1] < x->x[2] + x->x[3]);
	uint64_t mask = (1ULL<<k) - 1, y[2], hash;
	y[0] = bfc_hash_64(x->x[t<<1|0], mask);
	y[1] = bfc_hash_64(x->x[t<<1|1], mask);
	hash = (y[0] ^ y[1]) << k | ((y[0] + y[1]) & mask);
	ret = bfc_bf_insert(aux->bf, hash);
	if (ret == aux->opt.n_hashes) bfc_ch_insert(aux->ch, y);
}

static void worker(void *data, long k, int tid)
{
	bfc_aux_t *aux = (bfc_aux_t*)data;
	bseq1_t *s = &aux->seqs[k];
	const bfc_opt_t *o = &aux->opt;
	int i, l;
	bfc_kmer_t x = bfc_kmer_null;
	for (i = l = 0; i < s->l_seq; ++i) {
		int c = seq_nt6_table[(uint8_t)s->seq[i]] - 1;
		if (c < 4) {
			x.x[0] = (x.x[0]<<1 | (c&1))  & aux->mask;
			x.x[1] = (x.x[1]<<1 | (c>>1)) & aux->mask;
			x.x[2] = x.x[2]>>1 | (1ULL^(c&1))<<(o->k-1);
			x.x[3] = x.x[3]>>1 | (1ULL^c>>1) <<(o->k-1);
			if (++l >= o->k) bfc_kmer_insert(aux, &x);
		} else l = 0, x = bfc_kmer_null;
	}
}

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

int main(int argc, char *argv[])
{
	kseq_t *seq;
	gzFile fp;
	bfc_aux_t aux;
	int i, c;
	double t_real;

	t_real = realtime();
	bfc_opt_init(&aux.opt);
	while ((c = getopt(argc, argv, "k:s:b:L:t:h:")) >= 0) {
		if (c == 'k') aux.opt.k = atoi(optarg);
		else if (c == 'b') aux.opt.n_shift = atoi(optarg);
		else if (c == 't') aux.opt.n_threads = atoi(optarg);
		else if (c == 'h') aux.opt.n_hashes = atoi(optarg);
		else if (c == 'L' || c == 's') {
			double x;
			char *p;
			x = strtod(optarg, &p);
			if (*p == 'G' || *p == 'g') x *= 1e9;
			else if (*p == 'M' || *p == 'm') x *= 1e6;
			else if (*p == 'K' || *p == 'k') x *= 1e3;
			if (c == 's') {
				bfc_opt_by_size(&aux.opt, (long)x + 1);
				fprintf(stderr, "[M::%s] set k to %d\n", __func__, aux.opt.k);
			} else if (c == 'L') aux.opt.chunk_size = (long)x + 1;
		}
	}
	aux.mask = (1ULL<<aux.opt.k) - 1;

	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bfc [options] <in.fq>\n\n");
		fprintf(stderr, "Options: -s FLOAT     approx genome size (k/m/g allowed; change -k and -b) [unset]\n");
		fprintf(stderr, "         -k INT       k-mer length [%d]\n", aux.opt.k);
		fprintf(stderr, "         -t INT       number of threads [%d]\n", aux.opt.n_threads);
		fprintf(stderr, "         -b INT       set Bloom Filter size to pow(2,INT) bits [%d]\n", aux.opt.n_shift);
		fprintf(stderr, "         -h INT       use INT hash functions for Bloom Filter [%d]\n", aux.opt.n_hashes);
		fprintf(stderr, "\n");
		return 1;
	}

	aux.bf = bfc_bf_init(aux.opt.n_shift, aux.opt.n_hashes);
	aux.ch = bfc_ch_init(aux.opt.k);

	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while ((aux.seqs = bseq_read(seq, aux.opt.chunk_size, &aux.n_seqs)) != 0) {
		double rt, ct;
		rt = realtime(); ct = cputime();
		if (aux.opt.n_threads == 1)
			for (i = 0; i < aux.n_seqs; ++i)
				worker(&aux, i, 0);
		else kt_for(aux.opt.n_threads, worker, &aux, aux.n_seqs);
		for (i = 0; i < aux.n_seqs; ++i) {
			free(aux.seqs[i].seq); free(aux.seqs[i].qual);
		}
		free(aux.seqs);
		fprintf(stderr, "[M::%s] processed %d sequences in %.3f sec (%.1f%% CPU); # k-mers stored: %ld\n",
				__func__, aux.n_seqs, realtime() - rt, 100. * (cputime() - ct) / (realtime() - rt), (long)bfc_ch_count(aux.ch));
	}
	kseq_destroy(seq);
	gzclose(fp);

	bfc_ch_destroy(aux.ch);
	bfc_bf_destroy(aux.bf);

	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	return 0;
}
