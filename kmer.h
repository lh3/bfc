#ifndef BFC_KMER_H
#define BFC_KMER_H

typedef struct {
	uint64_t x[4];
} bfc_kmer_t;

static inline void bfc_kmer_append(int k, uint64_t x[4], int c)
{ // IMPORTANT: 0 <= c < 4
	uint64_t mask = (1ULL<<k) - 1;
	x[0] = (x[0]<<1 | (c&1))  & mask;
	x[1] = (x[1]<<1 | (c>>1)) & mask;
	x[2] = x[2]>>1 | (1ULL^(c&1))<<(k-1);
	x[3] = x[3]>>1 | (1ULL^c>>1) <<(k-1);
}

static inline void bfc_kmer_change(int k, uint64_t x[4], int d, int c) // d-bp from the 3'-end of k-mer; 0<=d<k
{ // IMPORTANT: 0 <= c < 4
	uint64_t t = ~(1ULL<<d);
	x[0] = (uint64_t) (c&1)<<d | (x[0]&t);
	x[1] = (uint64_t)(c>>1)<<d | (x[1]&t);
	t = ~(1ULL<<(k-1-d));
	x[2] = (uint64_t)(1^(c&1))<<(k-1-d) | (x[2]&t);
	x[3] = (uint64_t)(1^ c>>1)<<(k-1-d) | (x[3]&t);
}

#endif
