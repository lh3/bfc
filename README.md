BFC is a standalone high-performance tool for correcting sequencing errors from
Illumina sequencing data. It is specifically designed for high-coverage
whole-genome human data, though also performs well for small genomes.

The BFC algorithm is a variant of the classical spectrum alignment algorithm
introduced by Pevzner et al [PMID:11504945]. It uses an exhaustive search to
find a k-mer path through a read that minimizes a heuristic objective function
jointly considering penalties on correction, quality and k-mer support. This
algorithm was first implemented in my fermi assembler and then refined a few
times in fermi, fermi2 and now in bfc. In the k-mer counting phase, bfc uses a
blocked bloom filter to filter out most singleton k-mers and keeps the rest in
a hash table [PMID:21831268]. The use of bloom filter is how BFC is named
(though other correctors such as Lighter and Bless rely more on bloom filter
than BFC).

On 16 fast CPU cores, BFC is able to correct 1.6 billion 101bp reads in 6.5
wall clock hours, several times faster than other tools that can do the job
in 128GB RAM. It also appears to be more accurate based on preliminary
analyses. In particular, it rarely overcorrects (making a read worse than
the orginal) in comparison to others.

BFC can be invoked as:
```sh
bfc -s 3g -t16 reads.fq.gz | gzip -1 > corrected.fq.gz
```
where option `-s` specifies the approximate size of the genome. It is also
possible to use one set of reads to correct another set:
```sh
bfc -s 3g readset1.fq.gz readset2.fq.gz | gzip -1 > corrected_readset2.fq.gz
```
BFC also offers an option to trim reads containing unique k-mers (don't switch
`-s` and `-k` as some options are ordered):
```sh
bfc -1 -s 3g -k51 -t16 corrected.fq.gz | gzip -1 > trimmed.fq.gz
```
This command line keeps k-mer occuring twice or more in a bloom filter (with
some false positives) and identifies the longest stretch in a read that has
hits in the bloom filter. K-mer trimming is about four times as fast as error
correction.
