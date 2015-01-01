#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bfc.h"
#include "bbf.h"
#include "htab.h"
#include "bseq.h"

/* A note on multi-threading

   The bloom filter is always the same regardless of how many threads in use.
   However, the k-mer inserted to the hash table may be slightly different.
   Suppose k-mers A and B are both singletons and that if A is inserted first,
   B is a false positive and gets inserted to the hash table. In the
   multi-threading mode, nonetheless, B may be inserted before A. In this case,
   B is not a false positive any more. This is not a bug. The k-mers put into
   the hash table depends on the order of input.
*/

#define CNT_BUF_SIZE 256

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
double cputime(void);
double realtime(void);

typedef struct { // cache to reduce locking
	uint64_t y[2];
	int is_high;
} insbuf_t;

typedef struct {
	const bfc_opt_t *opt;
	bseq_file_t *ks;
	bfc_bf_t *bf, *bf_high;
	bfc_ch_t *ch;
	int *n_buf;
	insbuf_t **buf;
} bfc_cntaux_t;

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	bfc_cntaux_t *aux;
} bfc_count_data_t;

static int bfc_kmer_bufclear(bfc_cntaux_t *aux, int forced, int tid)
{
	int i, k, r;
	if (aux->ch == 0) return 0;
	for (i = k = 0; i < aux->n_buf[tid]; ++i) {
		r = bfc_ch_insert(aux->ch, aux->buf[tid][i].y, aux->buf[tid][i].is_high, forced);
		if (r < 0) aux->buf[tid][k++] = aux->buf[tid][i];
	}
	aux->n_buf[tid] = k;
	return k;
}

static void bfc_kmer_insert(bfc_cntaux_t *aux, const bfc_kmer_t *x, int is_high, int tid)
{
	int k = aux->opt->k, ret;
	uint64_t y[2], hash;
	hash = bfc_kmer_hash(k, x->x, y);
	ret = bfc_bf_insert(aux->bf, hash);
	if (ret == aux->opt->n_hashes) {
		if (aux->ch && bfc_ch_insert(aux->ch, y, is_high, 0) < 0) { // counting with a hash table
			insbuf_t *p;
			if (bfc_kmer_bufclear(aux, 0, tid) == CNT_BUF_SIZE)
				bfc_kmer_bufclear(aux, 1, tid);
			p = &aux->buf[tid][aux->n_buf[tid]++];
			p->y[0] = y[0], p->y[1] = y[1], p->is_high = is_high;
		} else if (aux->bf_high) // keep high-occurrence k-mers
			bfc_bf_insert(aux->bf_high, hash);
	}
}

static void worker_count(void *_data, long k, int tid)
{
	bfc_count_data_t *data = (bfc_count_data_t*)_data;
	bfc_cntaux_t *aux = data->aux;
	bseq1_t *s = &data->seqs[k];
	const bfc_opt_t *o = aux->opt;
	int i, l;
	bfc_kmer_t x = bfc_kmer_null;
	uint64_t qmer = 0, mask = (1ULL<<o->k) - 1;
	for (i = l = 0; i < s->l_seq; ++i) {
		int c = seq_nt6_table[(uint8_t)s->seq[i]] - 1;
		if (c < 4) {
			bfc_kmer_append(o->k, x.x, c);
			if (++l >= o->k) {
				qmer = (qmer<<1 | (s->qual == 0 || s->qual[i] - 33 >= o->q)) & mask;
				bfc_kmer_insert(aux, &x, (qmer == mask), tid);
			}
		} else l = 0, qmer = 0, x = bfc_kmer_null;
	}
}

static void *bfc_count_cb(void *shared, int step, void *_data)
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
		for (i = 0; i < aux->opt->n_threads; ++i)
			bfc_kmer_bufclear(aux, 1, i);
		for (i = 0; i < data->n_seqs; ++i) {
			bseq1_t *s = &data->seqs[i];
			free(s->seq); free(s->qual); free(s->name);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

void *bfc_count(const char *fn, const bfc_opt_t *opt)
{
	bfc_cntaux_t caux;
	void *ret;
	int i;

	memset(&caux, 0, sizeof(bfc_cntaux_t));
	caux.opt = opt;
	caux.bf = bfc_bf_init(opt->n_shift, opt->n_hashes);
	if (!opt->filter_mode) {
		caux.ch = bfc_ch_init(opt->k);
		caux.n_buf = calloc(opt->n_threads, sizeof(int));
		caux.buf = calloc(opt->n_threads, sizeof(void*));
		for (i = 0; i < opt->n_threads; ++i)
			caux.buf[i] = malloc(CNT_BUF_SIZE * sizeof(insbuf_t));
		caux.ks = bseq_open(fn);
		kt_pipeline(opt->no_mt_io? 1 : 2, bfc_count_cb, &caux, 2);
		bseq_close(caux.ks);
		for (i = 0; i < opt->n_threads; ++i) free(caux.buf[i]);
		free(caux.buf); free(caux.n_buf);
		ret = (void*)caux.ch;
	} else {
		caux.bf_high = bfc_bf_init(opt->n_shift, opt->n_hashes);
		caux.ks = bseq_open(fn);
		kt_pipeline(opt->no_mt_io? 1 : 2, bfc_count_cb, &caux, 2);
		bseq_close(caux.ks);
		ret = (void*)caux.bf_high;
	}
	bfc_bf_destroy(caux.bf);
	return ret;
}
