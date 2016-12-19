CFLAGS=-W -Wall -finline-functions -fPIC -std=gnu99 -Wno-unused-result

CLIB=-lpthread -lz -lm

OS := $(shell uname)
ifeq ($(OS),  Darwin)
  CFLAGS += -Wno-unused-function
else
	CLIB += -lrt
endif

INCLUDE = include

########### program ##########

# detect :
# 	@echo "$$CFLAGS" $(CFLAGS)

PROG = bin/biscuit
# PROG = bin/hemifinder bin/correct_bsstrand bin/get_unmapped bin/sample_trinuc

release : CFLAGS += -O3
release : $(PROG)

debug : CFLAGS += -g
debug : $(PROG)

######### libraries ###########

LHTSLIB_DIR = lib/htslib
LHTSLIB_INCLUDE = lib/htslib/htslib
LHTSLIB = $(LHTSLIB_DIR)/libhts.a
$(LHTSLIB) :
	make -C $(LHTSLIB_DIR) libhts.a

LKLIB_DIR = lib/klib
LKLIB = $(LKLIB_DIR)/klib.a
$(LKLIB) :
	make -C $(LKLIB_DIR) klib.a

LUTILS_DIR = lib/utils
LUTILS = $(LUTILS_DIR)/libutils.a
$(LUTILS):
	make -C $(LUTILS_DIR) libutils.a

LIBS=lib/aln/libaln.a src/pileup.o src/markdup.o src/ndr.o src/vcf2bed.o src/epiread.o $(LUTILS) $(LKLIB) $(LHTSLIB)
bin/biscuit: $(LIBS) src/main.o
	mkdir -p bin
	gcc $(CFLAGS) src/main.o -o $@ -I$(INCLUDE)/aln -I$(INCLUDE)/klib $(LIBS) $(CLIB)
clean_biscuit:
	rm -f bin/biscuit

####### subcommands #######

src/main.o: src/main.c
	gcc -c $(CFLAGS) src/main.c -o $@ -I$(LUTILS_DIR) -I$(LKLIB_DIR)
clean_main:
	rm -f src/main.o

LALND = lib/aln
LALNOBJ = $(LALND)/bntseq.o $(LALND)/bwamem.o $(LALND)/bwashm.o $(LALND)/bwt_gen.o $(LALND)/bwtsw2_chain.o $(LALND)/bwtsw2_pair.o $(LALND)/malloc_wrap.o $(LALND)/bwamem_extra.o $(LALND)/bwt.o $(LALND)/bwtindex.o $(LALND)/bwtsw2_core.o $(LALND)/fastmap.o  $(LALND)/QSufSort.o $(LALND)/bwa.o $(LALND)/bwamem_pair.o $(LALND)/bwtgap.o $(LALND)/bwtsw2_aux.o $(LALND)/bwtsw2_main.o $(LALND)/is.o $(LALND)/utils.o $(LALND)/ksw.o
lib/aln/libaln.a: $(LALNOBJ)
	ar -csru $@ $(LALNOBJ)
$(LALND)/%.o: $(LALND)/%.c
	gcc -c $(CFLAGS) -I$(LUTILS_DIR) -I$(INCLUDE)/klib $< -o $@
clean_aln:
	rm -f $(LALND)/*.o lib/aln/libaln.a

src/pileup.o: src/pileup.c
	gcc -c $(CFLAGS) -I$(LHTSLIB_INCLUDE) -I$(LUTILS_DIR) $< -o $@

src/markdup.o: src/markdup.c
	gcc -c $(CFLAGS) -I$(LHTSLIB_INCLUDE) -I$(LUTILS_DIR) $< -o $@

src/ndr.o: src/ndr.c
	gcc -c $(CFLAGS) -I$(LUTILS_DIR) -I$(LKLIB_DIR) $< -o $@

src/vcf2bed.o: src/vcf2bed.c
	gcc -c $(CFLAGS) -I$(LUTILS_DIR) -I$(LKLIB_DIR) $< -o $@

src/epiread.o: src/epiread.c
	gcc -c $(CFLAGS) -I$(LHTSLIB_INCLUDE) -I$(LUTILS_DIR) $< -o $@

####### general #######

.c.o :
	gcc -c $(CFLAGS) $< -o $@

####### clean #######

.PHONY: clean
clean :
	rm -f src/*.o

purge : clean
	make -C $(KLIB_DIR) clean
	make -C $(LHTSLIB) clean
	make -C $(LUTILS_DIR) clean

####### archived #######

.PHONY: correct_bsstrand
correct_bsstrand : bin/correct_bsstrand
bin/correct_bsstrand: $(LSAM0119)
	gcc $(CFLAGS) -o $@ -I$(INCLUDE) -I$(LSAM0119D) src/correct_bsstrand/correct_bsstrand.c $(LSAM0119) -lz -lpthread
clean_correct_bsstrand:
	rm -f bin/correct_bsstrand

# get trinuc spectrum
.PHONY: sample_trinuc
sample_trinuc : bin/sample_trinuc
bin/sample_trinuc: $(LSAM0119) src/sample_trinuc/sample_trinuc.c
	gcc $(CFLAGS) -o $@ -I$(LSAM0119D) -I$(INCLUDE) src/sample_trinuc/sample_trinuc.c -lpthread $(LSAM0119) $(LUTILS) -lz
clean_sample_trinuc:
	rm -f bin/sample_trinuc

# find hemi methylation
.PHONY: hemifinder
hemifinder : bin/hemifinder
bin/hemifinder: $(LSAM0119) $(LUTILS) src/hemifinder/hemifinder.c
	gcc $(CFLAGS) -o $@ -I$(LSAM0119D) -I$(INCLUDE) src/hemifinder/hemifinder.c $(LSAM0119) $(LUTILS) -lpthread  -lz
clean_hemifinder:
	rm -f bin/hemifinder
