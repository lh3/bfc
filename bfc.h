#ifndef BFC_H
#define BFC_H

#include "bbf.h"
#include "htab.h"
#include "bseq.h"

#define BFC_MAX_KMER     63
#define BFC_MAX_BF_SHIFT 37

#define BFC_MAX_PATHS 4
#define BFC_EC_HIST 5
#define BFC_EC_HIST_HIGH 2

typedef struct {
	int chunk_size;
	int n_threads, no_mt_io;
	int q, k;

	int filter_mode, refine_ec, no_qual;
	float min_frac;

	int l_pre, bf_shift, n_hashes;

	int discard;
	int max_end_ext;
	int win_multi_ec; // no 2 high-qual corrections or 4 corrections in a window of this size
	int min_cov; // a k-mer is considered solid if the count is no less than this

	// these ec options cannot be changed on the command line
	int w_ec, w_ec_high, w_absent, w_absent_high;
	int max_path_diff, max_heap;
} bfc_opt_t;

extern int bfc_verbose;
extern double bfc_real_time;
extern bfc_kmer_t bfc_kmer_null;

void *bfc_count(const char *fn, const bfc_opt_t *opt);
void bfc_correct(const char *fn, const bfc_opt_t *opt, const void *ptr);

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);
double cputime(void);
double realtime(void);

#endif
