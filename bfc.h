#ifndef BFC_H
#define BFC_H

#include "bbf.h"
#include "kmer.h"
#include "bseq.h"

#define BFC_MAX_KMER 63

#define BFC_MAX_PATHS 8
#define BFC_EC_HIST 5
#define BFC_EC_HIST_HIGH 2

typedef struct {
	int k, bf_shift; // modified based on the kmc counts

	int n_hashes; // number of hash functions for bloom filter
	int min_cov;

	int chunk_size;
	int n_threads, no_mt_io;

	int discard, no_qual;

	int q;
	int win_multi_ec;

	int kmer_trim;
	float trim_thres;

	// these ec options cannot be changed on the command line
	int max_end_ext;
	int w_ec, w_ec_high, w_absent;
	int max_path_diff, max_heap;
} bfc_opt_t;

extern int bfc_verbose;
extern double bfc_real_time;
extern bfc_kmer_t bfc_kmer_null;

bfc_bf_t *bfc_count(const char *fn, bfc_opt_t *opt);
void bfc_correct(const char *fn, const bfc_opt_t *opt, const bfc_bf_t *ptr);

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
double cputime(void);
double realtime(void);

#endif
