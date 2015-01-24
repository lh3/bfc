## Data

Raw reads were downloaded from [BaseSpace][basespace] under the sample
"NA12878-L7" of project "HiSeq X Ten: TruSeq Nano (4 replicates of NA12878)".
There are 444,764,118 pairs of reads, averaged 150bp or so in length. The raw
data are kept in two gzip'd files, one file for each read end. As some programs
work with one input file only or do not recognize gzip'd input files, we
may decompress the raw data or create one interleaved FASTQ.

## Hardware and Software

Fiona, BBMap and Bloocoo provide precompiled binaries. We compiled the rest of
tools on a CentOS5 virtual machine. The executables are [available
here][biobin].  We ran all the tools on a CentOS6 machine with 20 cores of
Intel E5-2660 CPUs at 2.2GHz and 128GB RAM, and measured timing and peak memory
with GNU time.

## Command Lines


```sh
# BBMap-34.38 (using ecclimit=14 made little difference)
ecc.sh -Xmx32g in=read1.fastq.gz in2=read2.fastq.gz out=ec.fq tmpdir=tmp threads=16 ecc=t aec=t k=31

# BFC-r155
bash -c "bfc -s 3g -k55 -t 16 <(seqtk mergepe read1.fq.gz read2.fq.gz) <(seqtk mergepe read1.fq.gz read2.fq.gz) | gzip -1 > ec.fq.gz"

# BFC-kmc
echo -e "read1.fq.gz\nread2.fq.gz" > list.txt
kmc -k55 -m24 @list.txt reads.k55 tmp
seqtk mergepe read1.fq.gz read2.fq.gz | bfc-kmc -t16 reads.k55 | gzip -1 > ec.fq.gz

# BLESS-v0p23 (much faster than v0p17 or earlier)
bless -read1 read1.fq -read2 read2.fq -kmerlength 55 -prefix out -smpthread 16 -max_mem 24 -notrim

# Bloocoo-1.0.4
Bloocoo -nb-cores 16 -file read12.fq -kmer-size 31

# Fermi2-r175; ropebwt2-r187
seqtk mergepe read1.fq.gz read2.fq.gz | ropebwt2 -drq20 -x31 > index.fmd
seqtk mergepe read1.fq.gz read2.fq.gz | fermi2 correct -t 16 -k 29 index.fmd /dev/stdin | gzip -1 > ec.fq.gz

# Fiona-0.2.0 (killed due to large memory footprint)
fiona -g 3000000000 --sequencing-technology illumina --no-final-trim-ns -nt 16 read12.fq ec.fq

# Lighter-20140123
lighter -K 31 3000000000 -r read1.fq.gz -r read2.fq.gz -t 16

# QuorUM-1.0.0
quorum -s 12000000000 -t 16 -p quorum-k31 -k 31 read1.fq read2.fq

# SGA-0.9.13
sga preprocess -p 1 read1.fq.gz read2.fq.gz | gzip -1 > out.pe.fq.gz
sga index -a ropebwt -t 16 --no-reverse out.pe.fq.gz
sga correct -t 16 -k 55 --learn out.pe.fq.gz

# Trowel-0.1.4.1 (killed due to large memory footprint)
echo "read1.fq read2.fq" > list.txt
trowel -t 16 -f list.txt -k 31 -ntr
```

## Comments

1. Early versions of BLESS were single-threaded and slow ([Molnar and Ilie,
   2014][review]), and did not work with unequal read lengths. More recent
   versions are much better.

2. In addition to tools shown in the table, we have also tried
   [Trowel-0.1.4.1][trowel], [Fiona-0.2.0][fiona] and
   [AllPathsLG-51828][allpath]. They were taking over 110GB RAM and got
   manually killed. [Coral][coral], [HiTEC][hitec] and [SHREC][shrec] are
   unable to work with data at this scale according to a few publications.
   [RACER][racer] requires more RAM than the total read bases (>120GB)
   [according to][review] the developers. These tools are not tested.

[basespace]: https://basespace.illumina.com/datacentral
[biobin]: https://sourceforge.net/projects/biobin/
[review]: http://bib.oxfordjournals.org/content/early/2014/09/01/bib.bbu029
[trowel]: https://sourceforge.net/projects/trowel-ec/
[fiona]: http://www.seqan.de/projects/fiona/
[allpath]: http://www.broadinstitute.org/software/allpaths-lg/blog/
[racer]: http://www.csd.uwo.ca/~ilie/RACER/
[time]: https://ftp.gnu.org/gnu/time/
[coral]: http://www.cs.helsinki.fi/u/lmsalmel/coral/
[hitec]: http://www.csd.uwo.ca/~ilie/HiTEC/
[shrec]: https://sourceforge.net/projects/shrec-ec/
