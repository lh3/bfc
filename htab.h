#ifndef BFC_HTAB_H
#define BFC_HTAB_H

#include <stdint.h>
#include "kmer.h"

struct bfc_ch_s;
typedef struct bfc_ch_s bfc_ch_t;

bfc_ch_t *bfc_ch_init(int k);
void bfc_ch_destroy(bfc_ch_t *ch);
int bfc_ch_insert(bfc_ch_t *ch, uint64_t x[2], int is_high, int forced);
int bfc_ch_get(const bfc_ch_t *ch, const uint64_t x[2]);
uint64_t bfc_ch_count(const bfc_ch_t *ch);
int bfc_ch_hist(const bfc_ch_t *ch, uint64_t cnt[256], uint64_t high[64]);
int bfc_ch_dump(const bfc_ch_t *ch, const char *fn);
bfc_ch_t *bfc_ch_restore(const char *fn);

struct kh_kc_s;
typedef struct kh_kc_s bfc_kc_t;

bfc_kc_t *bfc_kc_init(void);
void bfc_kc_destroy(bfc_kc_t *kc);
void bfc_kc_clear(bfc_kc_t *kc);
int bfc_kc_get(const bfc_ch_t *ch, bfc_kc_t *kc, const bfc_kmer_t *z);

// Thomas Wang's integer hash functions. See <https://gist.github.com/lh3/59882d6b96166dfc3d8d> for a snapshot.
static inline uint64_t bfc_hash_64(uint64_t key, uint64_t mask)
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

#endif
