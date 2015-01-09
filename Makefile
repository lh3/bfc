CC=			gcc
CFLAGS=		-g -Wall -O2 -Wno-unused-function #-fno-inline-functions -fno-inline-functions-called-once
CPPFLAGS=
INCLUDES=	
OBJS=		kthread.o utils.o bseq.o bbf.o htab.o count.o correct.o bfc.o
PROG=		bfc hash2cnt
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bfc:$(OBJS)
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

hash2cnt:hash2cnt.o
		$(CC) $(CFLAGS) $< -o $@

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bbf.o: bbf.h
bfc.o: bfc.h bbf.h htab.h kmer.h bseq.h
bseq.o: bseq.h kseq.h
correct.o: bfc.h bbf.h htab.h kmer.h bseq.h kvec.h ksort.h
count.o: bfc.h bbf.h htab.h kmer.h bseq.h
hash2cnt.o: kmer.h
htab.o: htab.h kmer.h khash.h
