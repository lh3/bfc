## Introduction

BFC is a standalone high-performance tool for correcting sequencing errors from
Illumina sequencing data. It is specifically designed for high-coverage
whole-genome human data, though also performs well for small genomes.

The BFC algorithm is a variant of the classical spectrum alignment algorithm
introduced by [Pevzner et al (2001)][Euler]. It uses an exhaustive search to
find a k-mer path through a read that minimizes a heuristic objective function
jointly considering penalties on correction, quality and k-mer support. This
algorithm was first implemented in my fermi assembler and then refined a few
times in fermi, fermi2 and now in BFC. In the k-mer counting phase, BFC uses a
blocked bloom filter to filter out most singleton k-mers and keeps the rest in a
hash table ([Melsted and Pritchard, 2011][bfcounter]). The use of bloom filter
is how BFC is named, though other correctors such as [Lighter][lighter] and
[Bless][bless] actually rely more on bloom filter than BFC.

## Usage

BFC can be invoked as:
```sh
bfc -s 3g -t16 reads.fq.gz | gzip -1 > corrected.fq.gz
```
where option `-s` specifies the approximate size of the genome. It is possible
to use one set of reads to correct another set:
```sh
bfc -s 3g -t16 readset1.fq.gz readset2.fq.gz | gzip -1 > corrected_readset2.fq.gz
```
and to process data from Unix pipes ("<(command)" is bash specific):
```sh
bash -c "bfc -s 3g -t16 <(bzip2 -dc reads.fq.bz2) <(bzip2 -dc reads.fq.bz2) | gzip -1 > out.fq.gz"
```
BFC also offers an option to trim reads containing singleton k-mers (don't switch
`-s` and `-k` as some options are ordered):
```sh
bfc -1 -s 3g -k51 -t16 corrected.fq.gz | gzip -1 > trimmed.fq.gz
```
This command line keeps k-mer occurring twice or more in a bloom filter (with
some false positives) and identifies the longest stretch in a read that has
hits in the bloom filter. K-mer trimming is about four times as fast as error
correction.

## BFC-KMC

An alternative implementation of the algorithm is available at the [kmc
branch][kmc-branch] of this repository. It uses [KMC2][kmc] for k-mer counting
and keeps high-occurrence k-mers in a bloom filter. BFC-KMC should be invoked as:
```sh
kmc -k55 reads.fq.gz prefix tmpdir
bfc-kmc -t16 prefix reads.fq.gz | gzip -1 > corrected.fq.gz
```
KMC2 source code and precompiled binaries are available at the [KMC
website][kmc].

[Euler]: http://www.ncbi.nlm.nih.gov/pubmed/11504945
[bfcounter]: http://www.ncbi.nlm.nih.gov/pubmed/21831268
[lighter]: https://github.com/mourisl/Lighter
[bless]: https://sourceforge.net/p/bless-ec/wiki/Home/
[kmc-branch]: https://github.com/lh3/bfc/tree/kmc
[kmc]: http://sun.aei.polsl.pl/kmc/index.html
