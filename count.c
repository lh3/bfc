#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bfc.h"
#include "bbf.h"
#include "bseq.h"
#include "bkmer.h"

typedef struct {
	const bfc_opt_t *opt;
	bkmer_file_t *bk;
	bfc_bf_t *bf;
} cnt_shared_t;

typedef struct {
	int n_kmers;
	bkmer1_t *kmers;
	cnt_shared_t *cs;
} cnt_step_t;

static void worker_count(void *_data, long k, int tid)
{
	cnt_step_t *data = (cnt_step_t*)_data;
	cnt_shared_t *cs = data->cs;
	char *kmer = data->kmers[k].kmer;
	uint64_t y[2], hash;
	int i, kk = cs->opt->k;
	bfc_kmer_t x;
	memset(x.x, 0, 32);
	for (i = 0; i < kk; ++i)
		bfc_kmer_append(kk, x.x, seq_nt6_table[(uint8_t)kmer[i]] - 1);
	hash = bfc_kmer_hash(kk, x.x, y);
	bfc_bf_insert(cs->bf, hash);
}

static void *bfc_count_cb(void *shared, int step, void *_data)
{
	cnt_shared_t *cs = (cnt_shared_t*)shared;
	if (step == 0) {
		cnt_step_t *ret;
		ret = calloc(1, sizeof(cnt_step_t));
		ret->kmers = bkmer_read(cs->bk, cs->opt->min_cov, cs->opt->chunk_size, &ret->n_kmers);
		ret->cs = cs;
		fprintf(stderr, "[M::%s] read %d k-mers\n", __func__, ret->n_kmers);
		if (ret->kmers) return ret;
		else free(ret);
	} else if (step == 1) {
		int i;
		double rt, eff;
		cnt_step_t *data = (cnt_step_t*)_data;
		kt_for(cs->opt->n_threads, worker_count, data, data->n_kmers);
		rt = realtime() - bfc_real_time;
		eff = 100. * cputime() / (rt + 1e-6);
		fprintf(stderr, "[M::%s @%.1f*%.1f%%] processed %d kmers\n", __func__, rt, eff, data->n_kmers);
		for (i = 0; i < data->n_kmers; ++i) free(data->kmers[i].kmer);
		free(data->kmers); free(data);
	}
	return 0;
}

bfc_bf_t *bfc_count(const char *fn, bfc_opt_t *opt)
{
	cnt_shared_t cs;
	uint64_t x;
	int i;
	memset(&cs, 0, sizeof(cnt_shared_t));
	cs.bk = bkmer_open(fn);
	opt->k = cs.bk->k;
	for (i = 0, x = 1; x < cs.bk->tot_kmers; ++i, x <<= 1);
	if (opt->bf_shift < i)
		opt->bf_shift = i + 4 < BFC_MAX_BF_SHIFT? i + 4 : BFC_MAX_BF_SHIFT;
	if (bfc_verbose >= 3)
		fprintf(stderr, "[M::%s] k=%d, tot_mers=%lld, bf_shift=%d\n", __func__, opt->k, (long long)cs.bk->tot_kmers, opt->bf_shift);
	cs.opt = opt;
	cs.bf = bfc_bf_init(opt->bf_shift, opt->n_hashes);
	kt_pipeline(opt->no_mt_io? 1 : 2, bfc_count_cb, &cs, 2);
	bkmer_close(cs.bk);
	return cs.bf;
}
