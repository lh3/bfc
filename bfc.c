#include <zlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

/**************
 * Base table *
 **************/

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

/******************
 * Hash functions *
 ******************/

// Thomas Wang's integer hash functions. See <https://gist.github.com/lh3/59882d6b96166dfc3d8d> for a snapshot.
uint64_t hash_64(uint64_t key, uint64_t mask)
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
uint64_t hash_64i(uint64_t key, uint64_t mask)
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

/*****************
 * Configuration *
 *****************/

typedef struct {
	int chunk_size;
	int n_threads;
	int k;
	int n_shift, n_hashes;
} bfc_opt_t;

void bfc_opt_init(bfc_opt_t *opt)
{
	memset(opt, 0, sizeof(bfc_opt_t));
	opt->chunk_size = 10000000;
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

typedef struct {
} bfc_hash_t;

/************************
 * Blocked Bloom Filter *
 ************************/

#define BFC_BLK_SHIFT  9 // 64 bytes, the size of a cache line
#define BFC_BLK_MASK   ((1<<(BFC_BLK_SHIFT)) - 1)

typedef struct {
	int n_shift, n_hashes;
	uint64_t *b;
} bfc_bf_t;

bfc_bf_t *bfc_bf_init(int n_shift, int n_hashes)
{
	bfc_bf_t *b;
	if (n_shift + BFC_BLK_SHIFT > 64 || n_shift < BFC_BLK_SHIFT) return 0;
	b = calloc(1, sizeof(bfc_bf_t));
	b->n_shift = n_shift;
	b->n_hashes = n_hashes;
	posix_memalign((void**)&b->b, 1<<(BFC_BLK_SHIFT-3), 1ULL<<(n_shift-3));
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
	uint64_t *p = &b->b[y<<(BFC_BLK_SHIFT-6)];
	int i, z = h1, cnt = 0;
	if (!(h2&1)) h2 = (h2 + 1) & BFC_BLK_MASK;
	for (i = z = 0; i < b->n_hashes; ++i) {
		uint64_t *q = &p[z>>6];
		uint64_t u = 1ULL<<(z&63);
		uint64_t v = __sync_fetch_and_or(q, u);
		if (v&u) ++cnt;
		z = (z + h2) & BFC_BLK_MASK;
	}
	return cnt;
}

int bfc_bf_test(bfc_bf_t *b, uint64_t hash)
{
	int x = b->n_shift - BFC_BLK_SHIFT;
	uint64_t y = hash & ((1ULL<<x) - 1);
	int h1 = hash >> x & BFC_BLK_MASK;
	int h2 = hash >> b->n_shift & BFC_BLK_MASK;
	uint64_t *p = &b->b[y<<(BFC_BLK_SHIFT-6)];
	int i, z = h1, cnt = 0;
	if (!(h2&1)) h2 = (h2 + 1) & BFC_BLK_MASK;
	for (i = z = 0; i < b->n_hashes; ++i) {
		if (p[z>>6] & 1ULL<<(z&63)) ++cnt;
		z = (z + h2) & BFC_BLK_MASK;
	}
	return cnt;
}

typedef struct {
	bfc_opt_t opt;
	uint64_t mask;
	int n_seqs;
	bseq1_t *seqs;
	bfc_bf_t *bf;
} bfc_aux_t;

void bfc_kmer_insert(bfc_aux_t *aux, const bfc_kmer_t *x)
{
	int k = aux->opt.k, ret;
	int t = !(x->x[0] + x->x[1] < x->x[2] + x->x[3]);
	uint64_t mask = (1ULL<<k) - 1, y[2], hash;
	y[0] = hash_64(x->x[t<<1|0], mask);
	y[1] = hash_64(x->x[t<<1|1], mask);
	hash = (y[0] ^ y[1]) << k | ((y[0] + y[1]) & mask);
	ret = bfc_bf_insert(aux->bf, hash);
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

	bfc_opt_init(&aux.opt);
	while ((c = getopt(argc, argv, "k:s:")) >= 0) {
		if (c == 'k') aux.opt.k = atoi(optarg);
		else if (c == 's') {
			long x;
			char *p;
			x = strtol(optarg, &p, 10);
			if (*p == 'G' || *p == 'g') x *= 1000000000;
			else if (*p == 'M' || *p == 'm') x *= 1000000;
			else if (*p == 'K' || *p == 'k') x *= 1000;
			bfc_opt_by_size(&aux.opt, x);
		}
	}
	aux.mask = (1ULL<<(aux.opt.k-1)) - 1;

	if (optind == argc) {
		fprintf(stderr, "Usage: bfc <in.fq>\n");
		return 1;
	}

	aux.bf = bfc_bf_init(aux.opt.n_shift, aux.opt.n_hashes);

	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while ((aux.seqs = bseq_read(seq, aux.opt.chunk_size, &aux.n_seqs)) != 0) {
		kt_for(aux.opt.n_threads, worker, &aux, aux.n_seqs);
		for (i = 0; i < aux.n_seqs; ++i) {
			free(aux.seqs[i].seq); free(aux.seqs[i].qual);
		}
		free(aux.seqs);
	}
	kseq_destroy(seq);
	gzclose(fp);

	bfc_bf_destroy(aux.bf);
	return 0;
}
