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
	int k, q;
	int n_shift, n_hashes;
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

int bfc_ch_get(const bfc_ch_t *ch, const bfc_kmer_t *z)
{
	int k = ch->k, ret = -1, t = !(z->x[0] + z->x[1] < z->x[2] + z->x[3]);
	uint64_t mask = (1ULL<<k) - 1, x[2], key;
	cnthash_t *h;
	khint_t itr;

	x[0] = bfc_hash_64(z->x[t<<1|0], mask);
	x[1] = bfc_hash_64(z->x[t<<1|1], mask);
	h = ch->h[x[0] & ((1ULL<<ch->l_pre) - 1)];
	key = (x[0] >> ch->l_pre | x[1] << (ch->k - ch->l_pre)) << 14 | 1;
	itr = kh_get(cnt, h, key);
	if (itr != kh_end(h)) {
		ret = kh_key(h, itr) & 0x3fff;
		if (!t) ret = (ret & 0xff) | (ret>>8&7)<<11 | (ret>>11&7)<<8;
	}
	return ret;
}

uint64_t bfc_ch_count(const bfc_ch_t *ch)
{
	int i;
	uint64_t cnt = 0;
	for (i = 0; i < 1<<ch->l_pre; ++i)
		cnt += kh_size(ch->h[i]);
	return cnt;
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

void _bfc_ch_print(const bfc_ch_t *ch)
{
	int i;
	for (i = 0; i < 1<<ch->l_pre; ++i) {
		cnthash_t *h = ch->h[i];
		khint_t k;
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k))
				printf("%d\t%d\t%d\n", (int)(kh_key(h, k)&0xff), (int)(kh_key(h,k)>>8&7), (int)(kh_key(h,k)>>11&7));
	}
}

void bfc_ch_print_kcov(const bfc_ch_t *ch, const char *seq)
{
	int len, i, l;
	uint64_t mask = (1ULL<<ch->k) - 1;
	len = strlen(seq);
	bfc_kmer_t x = bfc_kmer_null;
	for (i = l = 0; i < len; ++i) {
		int c = seq_nt6_table[(uint8_t)seq[i]] - 1;
		if (c < 4) {
			x.x[0] = (x.x[0]<<1 | (c&1))  & mask;
			x.x[1] = (x.x[1]<<1 | (c>>1)) & mask;
			x.x[2] = x.x[2]>>1 | (1ULL^(c&1))<<(ch->k-1);
			x.x[3] = x.x[3]>>1 | (1ULL^c>>1) <<(ch->k-1);
			if (++l >= ch->k) {
				int r;
				r = bfc_ch_get(ch, &x);
				if (r >= 0) fprintf(stderr, "%d\t%d\t%d\t%d\t%d\n", i - ch->k + 1, i + 1, r&0xff, r>>8&7, r>>11&7);
			}
		} else l = 0, x = bfc_kmer_null;
	}
}

/*******************************
 * Other routines for counting *
 *******************************/

static double bfc_real_time;

typedef struct {
	bfc_opt_t opt;
	uint64_t mask;
	kseq_t *ks;
	bfc_bf_t *bf;
	bfc_ch_t *ch;
} bfc_cnt_t;

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	bfc_cnt_t *aux;
} bfc_count_data_t;

void bfc_kmer_insert(bfc_cnt_t *aux, const bfc_kmer_t *x, int low_flag)
{
	int k = aux->opt.k, ret;
	int t = !(x->x[0] + x->x[1] < x->x[2] + x->x[3]);
	uint64_t mask = (1ULL<<k) - 1, y[2], hash;
	y[0] = bfc_hash_64(x->x[t<<1|0], mask);
	y[1] = bfc_hash_64(x->x[t<<1|1], mask);
	hash = (y[0] ^ y[1]) << k | ((y[0] + y[1]) & mask);
	ret = bfc_bf_insert(aux->bf, hash);
	if (t) low_flag = (low_flag&1)<<1 | (low_flag&2)>>1;
	if (ret == aux->opt.n_hashes)
		bfc_ch_insert(aux->ch, y, low_flag);
}

static void worker_count(void *_data, long k, int tid)
{
	bfc_count_data_t *data = (bfc_count_data_t*)_data;
	bfc_cnt_t *aux = data->aux;
	bseq1_t *s = &data->seqs[k];
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
	bfc_cnt_t *aux = (bfc_cnt_t*)shared;
	if (step == 0) {
		bfc_count_data_t *ret;
		ret = calloc(1, sizeof(bfc_count_data_t));
		ret->seqs = bseq_read(aux->ks, aux->opt.chunk_size, &ret->n_seqs);
		ret->aux = aux;
		fprintf(stderr, "[M::%s] read %d sequences\n", __func__, ret->n_seqs);
		if (ret->seqs) return ret;
		else free(ret);
	} else if (step == 1) {
		int i;
		bfc_count_data_t *data = (bfc_count_data_t*)_data;
		kt_for(aux->opt.n_threads, worker_count, data, data->n_seqs);
		fprintf(stderr, "[M::%s] processed %d sequences (CPU/real time: %.3f/%.3f secs; # distinct k-mers: %ld)\n",
				__func__, data->n_seqs, cputime(), realtime() - bfc_real_time, (long)bfc_ch_count(aux->ch));
		for (i = 0; i < data->n_seqs; ++i) {
			free(data->seqs[i].seq); free(data->seqs[i].qual);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

int main(int argc, char *argv[])
{
	gzFile fp;
	bfc_cnt_t aux;
	int i, c, no_mt_io = 0;
	char *in_hash = 0, *out_hash = 0, *str_kcov = 0;

	bfc_real_time = realtime();
	bfc_opt_init(&aux.opt);
	while ((c = getopt(argc, argv, "d:k:s:b:L:t:C:h:q:Jr:")) >= 0) {
		if (c == 'k') aux.opt.k = atoi(optarg);
		else if (c == 'C') str_kcov = optarg;
		else if (c == 'd') out_hash = optarg;
		else if (c == 'r') in_hash = optarg;
		else if (c == 'q') aux.opt.q = atoi(optarg);
		else if (c == 'b') aux.opt.n_shift = atoi(optarg);
		else if (c == 't') aux.opt.n_threads = atoi(optarg);
		else if (c == 'h') aux.opt.n_hashes = atoi(optarg);
		else if (c == 'J') no_mt_io = 1; // for debugging kt_pipeline()
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
		fprintf(stderr, "         -d FILE      dump hash table to FILE [null]\n");
		fprintf(stderr, "         -r FILE      restore hash table from FILE [null]\n");
		fprintf(stderr, "\n");
		return 1;
	}

	if (in_hash == 0) {
		aux.bf = bfc_bf_init(aux.opt.n_shift, aux.opt.n_hashes);
		aux.ch = bfc_ch_init(aux.opt.k);

		fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
		aux.ks = kseq_init(fp);
		kt_pipeline(no_mt_io? 1 : 2, bfc_count_cb, &aux, 2);
		kseq_destroy(aux.ks);
		gzclose(fp);

		bfc_bf_destroy(aux.bf);
		aux.bf = 0;
	} else aux.ch = bfc_ch_restore(in_hash);

	if (str_kcov) bfc_ch_print_kcov(aux.ch, str_kcov);
	if (out_hash) bfc_ch_dump(aux.ch, out_hash);
	bfc_ch_destroy(aux.ch);
	bfc_bf_destroy(aux.bf);

	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - bfc_real_time, cputime());
	return 0;
}
