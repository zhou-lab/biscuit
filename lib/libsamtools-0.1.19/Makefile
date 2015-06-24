# CFLAGS=-finline-functions -fPIC -g # no -W to suppress Warning from samtools
CFLAGS=-finline-functions -fPIC -O3 # -g no -W to suppress Warning from samtools

LOBJ = bam.o sam_header.o bgzf.o kstring.o bam_aux.o sam.o bam_import.o faidx.o razf.o bam_pileup.o bam_index.o bam_sort.o bam_rmdup.o bam_rmdupse.o knetfile.o

libsam.a : $(LOBJ)
	ar -csru $@ $(LOBJ)

samtools/%.o : samtools/%.c
	gcc -c $(CFLAGS) $< -o $@

test : 
	gcc bam.c main.c sam_header.c bgzf.c kstring.c bam_aux.c sam.c bam_import.c faidx.c razf.c bam_pileup.c -lpthread -lz -Wall -g -o test

clean :
	rm -fr gmon.out *.o a.out *.exe *.dSYM razip bgzip *~ *.a *.so.* *.so *.dylib
