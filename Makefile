CFLAGS=-W -Wall -finline-functions -fPIC -std=gnu99

OS := $(shell uname)
ifeq ($(OS),  Darwin)
  CFLAGS += -Wno-unused-function
endif

INCLUDE = include

LSAM0119D = lib/libsamtools-0.1.19
LSAM0119 = $(LSAM0119D)/libsam.a

LUTILSD = lib/utils
LUTILS = lib/utils/libutils.a

PROG = biscuit
# PROG = hemifinder correct_bsstrand get_unmapped sample_trinuc

release : CFLAGS += -O3
release : $(PROG)

debug : CFLAGS += -g
debug : $(PROG)

$(LSAM0119) :
	make -C $(LSAM0119D) libsam.a

.PHONY: klib
klib: lib/klib/klib.a
KLIBD = lib/klib
KLIBOBJ = $(KLIBD)/kstring.o $(KLIBD)/kopen.o $(KLIBD)/kthread.o $(KLIBD)/ksw.o
lib/klib/klib.a: $(KLIBOBJ)
	ar -csru $@ $(KLIBOBJ)
$(KLIBD)/%.o: $(KLIBD)/%.c
	gcc -c $(CFLAGS) -I$(INCLUDE)/klib $< -o $@
clean_klib:
	rm -f $(KLIBD)/*.o lib/klib/klib.a

.PHONY: biscuit
LIBS=lib/aln/libaln.a src/pileup.o src/somatic.o lib/klib/klib.a $(LSAM0119) $(LUTILS)
biscuit: bin/biscuit
bin/biscuit: src/main.c $(LIBS)
		$(CC) $(CFLAGS) src/main.c -o $@ -I$(INCLUDE)/aln -I$(INCLUDE)/klib $(LIBS) -lpthread -lz -lm -lrt

.PHONY: aln
aln: lib/aln/libaln.a
LALND = lib/aln
LALNOBJ = $(LALND)/bntseq.o $(LALND)/bwamem.o $(LALND)/bwashm.o $(LALND)/bwt_gen.o $(LALND)/bwtsw2_chain.o $(LALND)/bwtsw2_pair.o $(LALND)/malloc_wrap.o $(LALND)/bwamem_extra.o $(LALND)/bwt.o $(LALND)/bwtindex.o $(LALND)/bwtsw2_core.o $(LALND)/fastmap.o  $(LALND)/QSufSort.o $(LALND)/bwa.o $(LALND)/bwamem_pair.o $(LALND)/bwtgap.o $(LALND)/bwtsw2_aux.o $(LALND)/bwtsw2_main.o $(LALND)/is.o $(LALND)/utils.o
lib/aln/libaln.a: $(LALNOBJ)
	ar -csru $@ $(LALNOBJ)
$(LALND)/%.o: $(LALND)/%.c
	gcc -c $(CFLAGS) -I$(INCLUDE)/aln -I$(INCLUDE)/klib $< -o $@
clean_aln:
	rm -f $(LALND)/*.o lib/aln/libaln.a

.PHONY: utils
utils: $(LUTILS)
LUTILSOBJ = $(LUTILSD)/encode.o $(LUTILSD)/stats.o
$(LUTILS): $(LUTILSOBJ)
	ar -csru $@ $(LUTILSOBJ)
$(LUTILSD)/%.o: $(LUTILSD)/%.c
	gcc -c $(CFLAGS) -I$(INCLUDE) $< -o $@
clean_utils:
	rm -f $(LUTILSD)/*.o $(LUTILS)

libaln.a: $(ALNOBJS)
	$(AR) -csru $@ $(ALNOBJS)

.PHONY: pileup
pileup: src/pileup.o
src/pileup.o: src/pileup.c
	gcc -c $(CFLAGS) -o $@ -I$(LSAM0119D) -I$(INCLUDE) src/pileup.c
clean_pileup:
	rm -f src/pileup.o

.PHONY: somatic
somatic: src/somatic.o
somatic: src/somatic.o
src/somatic.o: src/somatic.c
	gcc -c $(CFLAGS) -o $@ -I$(LSAM0119D) -I$(INCLUDE) src/somatic.c
clean_somatic:
	rm -f src/somatic.o

.PHONY: correct_bsstrand
correct_bsstrand : bin/correct_bsstrand
bin/correct_bsstrand: $(LSAM0119)
	gcc $(CFLAGS) -o $@ -I$(INCLUDE) -I$(LSAM0119D) src/correct_bsstrand/correct_bsstrand.c $(LSAM0119) -lz -lpthread
clean_correct_bsstrand:
	rm -f bin/correct_bsstrand

.PHONY: get_unmapped
get_unmapped : bin/get_unmapped
bin/get_unmapped : $(LSAM0119)
	gcc $(CFLAGS) -o $@ -I$(LSAM0119D) src/get_unmapped/get_unmapped.c $(LSAM0119) -lz -lpthread
clean_get_unmapped:
	rm -f bin/get_unmapped

.PHONY: sample_trinuc
sample_trinuc : bin/sample_trinuc
bin/sample_trinuc: $(LSAM0119) src/sample_trinuc/sample_trinuc.c
	gcc $(CFLAGS) -o $@ -I$(LSAM0119D) -I$(INCLUDE) src/sample_trinuc/sample_trinuc.c -lpthread $(LSAM0119) $(LUTILS) -lz
clean_sample_trinuc:
	rm -f bin/sample_trinuc

.PHONY: hemifinder
hemifinder : bin/hemifinder
bin/hemifinder: $(LSAM0119) $(LUTILS) src/hemifinder/hemifinder.c
	gcc $(CFLAGS) -o $@ -I$(LSAM0119D) -I$(INCLUDE) src/hemifinder/hemifinder.c $(LSAM0119) $(LUTILS) -lpthread  -lz
clean_hemifinder:
	rm -f bin/hemifinder

.c.o :
	gcc -c $(CFLAGS) $< -o $@

.PHONY: clean_all
clean_all : clean_sample_trinuc clean_get_unmapped clean_correct_bsstrand clean_hemifinder clean_pileup_cytosine clean_klib
	make -C $(LSAM0119D) clean
