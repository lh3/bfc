#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include "bfc.h"

#define BFC_VERSION "r150"

int bfc_verbose = 3;
double bfc_real_time;
bfc_kmer_t bfc_kmer_null = {{0,0,0,0}};

void bfc_opt_init(bfc_opt_t *opt)
{
	memset(opt, 0, sizeof(bfc_opt_t));
	opt->chunk_size = 100000000;
	opt->n_threads = 1;
	opt->q = 20;
	opt->n_hashes = 4;

	opt->min_cov = 4;
	opt->win_multi_ec = 10;

	opt->max_end_ext = 5;
	opt->w_ec = 1;
	opt->w_ec_high = 7;
	opt->w_absent = 3;
	opt->max_path_diff = 15;
	opt->max_heap = 100;
}

static void usage(FILE *fp, bfc_opt_t *o)
{
	fprintf(fp, "Usage: bfc [options] <kmc.prefix>|<bf.dump> [in.fq]\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -t INT       number of threads [%d]\n", o->n_threads);
	fprintf(fp, "  -H INT       use INT hash functions for Bloom filter [%d]\n", o->n_hashes);
	fprintf(fp, "  -c INT       min k-mer coverage [%d]\n", o->min_cov);
	fprintf(fp, "  -w INT       no more than %d ec or %d highQ ec in INT-bp window [%d]\n", BFC_EC_HIST, BFC_EC_HIST_HIGH, o->win_multi_ec);
//	fprintf(fp, "  -D           discard uncorrectable reads\n");
	fprintf(fp, "  -Q           force FASTA output\n");
	fprintf(fp, "  -v           show version number\n");
	fprintf(fp, "  -h           show command line help\n");
}

int main(int argc, char *argv[])
{
	bfc_opt_t opt;
	bfc_bf_t *bf = 0;
	int i, c;

	bfc_real_time = realtime();
	bfc_opt_init(&opt);
	while ((c = getopt(argc, argv, "hvV:k:L:t:H:q:Jc:w:DQ")) >= 0) {
		if (c == 'q') opt.q = atoi(optarg);
		else if (c == 't') opt.n_threads = atoi(optarg);
		else if (c == 'H') opt.n_hashes = atoi(optarg);
		else if (c == 'c') opt.min_cov = atoi(optarg);
		else if (c == 'w') opt.win_multi_ec = atoi(optarg);
		else if (c == 'D') opt.discard = 1;
		else if (c == 'Q') opt.no_qual = 1;
		else if (c == 'J') opt.no_mt_io = 1; // for debugging kt_pipeline()
		else if (c == 'V') bfc_verbose = atoi(optarg);
		else if (c == 'h') {
			usage(stdout, &opt);
			return 0;
		} else if (c == 'v') {
			printf("%s\n", BFC_VERSION);
			return 0;
		} else if (c == 'L') {
			double x;
			char *p;
			x = strtod(optarg, &p);
			if (*p == 'G' || *p == 'g') x *= 1e9;
			else if (*p == 'M' || *p == 'm') x *= 1e6;
			else if (*p == 'K' || *p == 'k') x *= 1e3;
			opt.chunk_size = (long)x + 1;
		}
	}

	if (optind == argc) {
		usage(stderr, &opt);
		return 1;
	}

	if (bfc_bf_is_dump(argv[optind])) bf = bfc_bf_restore(argv[optind], &opt.k);
	else bf = bfc_count(argv[optind], &opt);
	if (optind + 1 < argc) bfc_correct(argv[optind+1], &opt, bf);
	else bfc_bf_dump(0, opt.k, bf);
	bfc_bf_destroy(bf);

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, BFC_VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - bfc_real_time, cputime());
	return 0;
}
