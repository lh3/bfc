#include "bkmer.h"
#include "kmc_file.h"

bkmer_file_t *bkmer_open(const char *fn)
{
	bkmer_file_t *r;
	CKMCFile *kmc;
	kmc = new CKMCFile();
	kmc->OpenForListing(fn);
	r = (bkmer_file_t*)calloc(1, sizeof(bkmer_file_t));
	r->kmc = (void*)kmc;
	r->k = kmc->KmerLength();
	r->tot_kmers = kmc->KmerCount();
	return r;
}

bkmer1_t *bkmer_read(bkmer_file_t *fp, int min_occ, int chunk_size, int *n_)
{
	CKMCFile *kmc = (CKMCFile*)fp->kmc;
	CKmerAPI kmer(fp->k);
	float count;
	int size = 0, m = 0, n = 0;
	bkmer1_t *kmers = 0;
	while (kmc->ReadNextKmer(kmer, count)) {
		bkmer1_t *s;
		int c = (int)(count + .499);
		if (c < min_occ) continue;
		if (n >= m) {
			m = m? m<<1 : 256;
			kmers = (bkmer1_t*)realloc(kmers, m * sizeof(bkmer1_t));
		}
		s = &kmers[n++];
		s->kmer = (char*)malloc(fp->k + 1);
		kmer.to_string(s->kmer);
		s->count = c;
		size += fp->k;
		if (size >= chunk_size) break;
	}
	*n_ = n;
	return kmers;
}

void bkmer_close(bkmer_file_t *fp)
{
	CKMCFile *kmc = (CKMCFile*)fp->kmc;
	delete kmc;
	free(fp);
}
