CFLAGS=		-g -Wall -O2 -Wc++-compat #-Wextra
OBJS=
PROG=		gen-fq
LIBS=		-lm -lz -lpthread

.PHONY:all clean depend
.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

gen-fq:gen-fq.o
		$(CC) -o $@ $< $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

gen-fq.o: ketopt.h kseq.h khash.h
