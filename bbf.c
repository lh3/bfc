#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bbf.h"

#define BBF_MAGIC "BBF\1"
#define BBF_READ_SIZE 0x1000000

bfc_bf_t *bfc_bf_init(int n_shift, int n_hashes)
{
	bfc_bf_t *b;
	void *ptr = 0;
	if (n_shift + BFC_BLK_SHIFT > 64 || n_shift < BFC_BLK_SHIFT) return 0;
	b = calloc(1, sizeof(bfc_bf_t));
	b->n_shift = n_shift;
	b->n_hashes = n_hashes;
	posix_memalign(&ptr, 1<<(BFC_BLK_SHIFT-3), 1ULL<<(n_shift-3));
	b->b = ptr;
	bzero(b->b, 1ULL<<(n_shift-3));
	return b;
}

void bfc_bf_destroy(bfc_bf_t *b)
{
	if (b == 0) return;
	free(b->b); free(b);
}

int bfc_bf_insert(bfc_bf_t *b, uint64_t hash)
{
	int x = b->n_shift - BFC_BLK_SHIFT;
	uint64_t y = hash & ((1ULL<<x) - 1);
	int h1 = hash >> x & BFC_BLK_MASK;
	int h2 = hash >> b->n_shift & BFC_BLK_MASK;
	uint8_t *p = &b->b[y<<(BFC_BLK_SHIFT-3)];
	int i, z = h1, cnt = 0;
	if ((h2&31) == 0) h2 = (h2 + 1) & BFC_BLK_MASK; // otherwise we may repeatedly use a few bits
	while (__sync_lock_test_and_set(p, 1)); // lock
	for (i = 0; i < b->n_hashes; z = (z + h2) & BFC_BLK_MASK) {
		uint8_t *q = &p[z>>3], u;
		if (p == q) continue; // don't use the first byte. It is a spin lock.
		u = 1<<(z&7);
		cnt += !!(*q & u);
		*q |= u;
		++i;
	}
	__sync_lock_release(p); // unlock
	return cnt;
}

int bfc_bf_get(const bfc_bf_t *b, uint64_t hash)
{
	int x = b->n_shift - BFC_BLK_SHIFT;
	uint64_t y = hash & ((1ULL<<x) - 1);
	int h1 = hash >> x & BFC_BLK_MASK;
	int h2 = hash >> b->n_shift & BFC_BLK_MASK;
	uint8_t *p = &b->b[y<<(BFC_BLK_SHIFT-3)];
	int i, z = h1, cnt = 0;
	if ((h2&31) == 0) h2 = (h2 + 1) & BFC_BLK_MASK; // otherwise we may repeatedly use a few bits
	for (i = 0; i < b->n_hashes; z = (z + h2) & BFC_BLK_MASK) {
		uint8_t *q = &p[z>>3];
		if (p == q) continue; // don't use the first byte. It is a spin lock.
		cnt += !!(*q & 1<<(z&7));
		++i;
	}
	return cnt;
}

void bfc_bf_dump(const char *fn, const bfc_bf_t *bf)
{
	FILE *fp;
	uint32_t x[2];
	uint64_t rest;
	uint8_t *p;
	fp = fn && strcmp(fn, "-")? fopen(fn, "wb") : stdout;
	if (fp == 0) return;
	x[0] = bf->n_shift, x[1] = bf->n_hashes;
	fwrite(BBF_MAGIC, 1, 4, fp);
	fwrite(x, 4, 2, fp);
	rest = 1ULL<<(bf->n_shift-3), p = bf->b;
	while (rest) {
		size_t ret;
		ret = fwrite(p, 1, BBF_READ_SIZE, fp);
		p += ret; rest -= ret;
	}
	fclose(fp);
}

int bfc_bf_is_dump(const char *fn)
{
	FILE *fp;
	char magic[4];
	if ((fp = fopen(fn, "rb")) == 0) return 0;
	fread(magic, 1, 4, fp);
	fclose(fp);
	return (strncmp(magic, BBF_MAGIC, 4) == 0);
}

bfc_bf_t *bfc_bf_restore(const char *fn)
{
	FILE *fp;
	uint32_t x[2];
	char magic[4];
	uint64_t rest;
	bfc_bf_t *bf;
	uint8_t *p;

	fp = fn && strcmp(fn, "-")? fopen(fn, "rb") : stdin;
	if (fp == 0) return 0;
	fread(magic, 1, 4, fp);
	if (strncmp(magic, BBF_MAGIC, 4)) {
		fclose(fp);
		return 0;
	}
	fread(x, 4, 2, fp);
	bf = bfc_bf_init(x[0], x[1]);
	rest = 1ULL<<(bf->n_shift-3), p = bf->b;
	while (rest) {
		size_t ret;
		ret = fread(p, 1, BBF_READ_SIZE, fp);
		p += ret; rest -= ret;
	}
	fclose(fp);
	return bf;
}
