#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include "kmer.h"

int main(int argc, char *argv[])
{
	int c, mode = 2, l_pre, k, i, j;
	uint32_t t[2];
	uint64_t mask;
	char buf[64];
	FILE *fp;
	while ((c = getopt(argc, argv, "m:")) >= 0) {
		if (c == 'm') mode = atoi(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: hash2cnt [options] <dump.hash>\n");
		return 1;
	}
	if ((fp = fopen(argv[optind], "rb")) == 0) return 1;
	fread(t, 4, 2, fp);
	k = t[0], l_pre = t[1];
	mask = (1ULL<<k) - 1;
	for (i = 0; i < 1<<l_pre; ++i) {
		uint64_t tmp;
		fread(t, 4, 2, fp);
		if (mode == 1)
			printf("%d\n", t[1]);
		for (j = 0; j < t[1]; ++j) {
			fread(&tmp, 8, 1, fp);
			if (mode == 2 || mode == 3) {
				uint64_t h[2], y[2];
				uint16_t cnt = tmp & 0x3fff;
				if (k <= 32) {
					uint64_t z = (uint64_t)i<<(k*2-l_pre) | tmp>>14;
					h[0] = z >> k;
					h[1] = z & mask;
				} else {
					h[0] = (uint64_t)i<<(k-l_pre) | tmp>>(14+k);
					h[1] = (tmp>>14) & mask;
				}
				bfc_kmer_hash_inv(k, h, y);
				if (mode == 2) {
					printf("%s\t%d\t%d\n", bfc_kmer_2str(k, y, buf), cnt&0xff, cnt>>8&0x3f);
				} else {
					fwrite(h, 8, 2, stdout);
					fwrite(&cnt, 2, 1, stdout);
				}
			}
		}
	}
	fclose(fp);
	return 0;
}
