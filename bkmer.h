#ifndef BKMER_H
#define BKMER_H

#include <stdint.h>

typedef struct {
	void *kmc;
	uint64_t k:8, tot_kmers:56;
} bkmer_file_t;

typedef struct {
	char *kmer;
	int count;
} bkmer1_t;

#ifdef __cplusplus
extern "C" {
#endif

bkmer_file_t *bkmer_open(const char *fn);
bkmer1_t *bkmer_read(bkmer_file_t *fp, int min_occ, int chunk_size, int *n_);
void bkmer_close(bkmer_file_t *fp);

#ifdef __cplusplus
}
#endif

#endif
