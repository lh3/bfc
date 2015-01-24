#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include "kmer.h"

int main(int argc, char *argv[])
{
	int c, subcnt_only = 0, hist_only = 0, l_pre, k, i, j;
	int min_cnt = 0, diff_thres = 0;
	uint32_t t[2];
	uint64_t mask, hist_all[256], hist_high[64];
	char buf[64];
	FILE *fp;

	while ((c = getopt(argc, argv, "shm:d:")) >= 0) {
		if (c == 's') subcnt_only = 1;
		else if (c == 'h') hist_only = 1;
		else if (c == 'm') min_cnt = atoi(optarg);
		else if (c == 'd') diff_thres = atoi(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: hash2cnt [options] <dump.hash>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -s       only show # elements in each sub- hash table\n");
		fprintf(stderr, "  -h       only show k-mer histogram\n");
		fprintf(stderr, "  -m INT   occ >= INT [0]\n");
		fprintf(stderr, "  -d INT   occ - occHigh >= INT [0]\n");
		return 1;
	}
	memset(hist_all, 0, 256*8);
	memset(hist_high,0, 64*8);

	if ((fp = fopen(argv[optind], "rb")) == 0) return 1;
	fread(t, 4, 2, fp);
	k = t[0], l_pre = t[1];
	if (k > 37) {
		fprintf(stderr, "ERROR: hash2cnt does not work for k>37\n");
		fclose(fp);
		return 1;
	}
	mask = (1ULL<<k) - 1;
	for (i = 0; i < 1<<l_pre; ++i) {
		uint64_t tmp;
		fread(t, 4, 2, fp);
		if (subcnt_only) printf("%d\n", t[1]);
		for (j = 0; j < t[1]; ++j) {
			int high, all, diff;
			fread(&tmp, 8, 1, fp);
			high = (tmp>>8&0x3f), all = tmp&0xff;
			diff = (all < 0x3f? all : 0x3f) - high;
			++hist_all[all]; ++hist_high[high];
			if (!subcnt_only && !hist_only && all >= min_cnt && diff >= diff_thres) {
				uint64_t h[2], y[2];
				if (k <= 32) {
					uint64_t z = (uint64_t)i<<(k*2-l_pre) | tmp>>14;
					h[0] = z >> k;
					h[1] = z & mask;
				} else {
					h[0] = (uint64_t)i<<(k-l_pre) | tmp>>(14+k);
					h[1] = (tmp>>14) & mask;
				}
				bfc_kmer_hash_inv(k, h, y);
				printf("%s\t%d\t%d\n", bfc_kmer_2str(k, y, buf), all, high);
			}
		}
	}
	fclose(fp);

	if (hist_only) {
		for (i = 0; i <256; ++i)
			if (i >= 64) printf("%d\t%llu\n", i, (unsigned long long)hist_all[i]);
			else printf("%d\t%llu\t%llu\n", i, (unsigned long long)hist_all[i], (unsigned long long)hist_high[i]);
	}
	return 0;
}
