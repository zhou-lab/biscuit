CC = gcc
CFLAGS = -W -Wall -finline-functions -fPIC -std=gnu99 -Wno-unused-result -O3
CLIB = -lncurses -lpthread -lz -lm
CF_OPTIMIZE = 1

OS := $(shell uname)
ifeq ($(OS),  Darwin)
	CFLAGS += -Wno-unused-function
else
	CLIB += -lrt -ltinfo
endif

INCLUDE = include

########################
### different modes ####
########################

# detect :
# 	@echo "$$CFLAGS" $(CFLAGS)

PROG = biscuit

.PHONY : setdebug debug build

# ifeq (1, $(CF_NO_OPTIMIZE))
#   CFLAGS += -g
# else
#   CFLAGS += -O3
# endif

build: exportcf $(PROG)

debug: CF_OPTIMIZE := 0
debug: CFLAGS += -g
debug: CFLAGS := $(filter-out -O3,$(CFLAGS))
debug: build

exportcf:
	$(eval export CF_OPTIMIZE)

#####################
##### libraries #####
#####################

LHTSLIB_DIR = lib/htslib
LHTSLIB_INCLUDE = lib/htslib/htslib
LHTSLIB = $(LHTSLIB_DIR)/libhts.a
$(LHTSLIB) :
	make -C $(LHTSLIB_DIR) libhts.a

LKLIB_DIR = lib/klib
LKLIB = $(LKLIB_DIR)/klib2.a
$(LKLIB) :
	make -C $(LKLIB_DIR) klib2.a

LUTILS_DIR = lib/utils
LUTILS = $(LUTILS_DIR)/libutils.a
$(LUTILS):
	make -C $(LUTILS_DIR) libutils.a

LSGSL_DIR = lib/sgsl
LSGSL = $(LSGSL_DIR)/libgsl.a
$(LSGSL):
	make -C $(LSGSL_DIR) libgsl.a

###########################
###### main program #######
###########################

BISCUITSRCS := $(wildcard src/*.c)
BISCUITLIBS := $(BISCUITSRCS:.c=.o)

LIBS=lib/aln/libaln.a src/pileup.o src/markdup.o src/ndr.o src/vcf2bed.o src/epiread.o src/asm_pairwise.o src/tview.o src/bsstrand.o src/cinread.o src/mergecg.o src/bsconv.o src/bamfilter.o $(LUTILS) $(LKLIB) $(LHTSLIB) $(LSGSL)
biscuit: $(LIBS) src/main.o
	gcc $(CFLAGS) src/main.o -o $@ -I$(INCLUDE)/aln -I$(INCLUDE)/klib $(LIBS) $(CLIB)

clean_biscuit:
	rm -f biscuit

###################
### subcommands ###
###################

src/main.o: src/main.c
	gcc -c $(CFLAGS) src/main.c -o $@ -I$(LUTILS_DIR) -I$(LKLIB_DIR)

LALND = lib/aln
LALNOBJ=$(patsubst %.c,%.o,$(wildcard $(LALND)/*.c))
# LALNOBJ = $(LALND)/bntseq.o $(LALND)/bwamem.o $(LALND)/bwashm.o $(LALND)/bwt_gen.o $(LALND)/bwtsw2_chain.o $(LALND)/bwtsw2_pair.o $(LALND)/malloc_wrap.o $(LALND)/bwamem_extra.o $(LALND)/bwt.o $(LALND)/bwtindex.o $(LALND)/bwtsw2_core.o $(LALND)/align.o  $(LALND)/QSufSort.o $(LALND)/bwa.o $(LALND)/bwamem_pair.o $(LALND)/bwtgap.o $(LALND)/bwtsw2_aux.o $(LALND)/bwtsw2_main.o $(LALND)/is.o $(LALND)/utils.o $(LALND)/ksw.o $(LALND)/mem_pair.o
lib/aln/libaln.a: $(LALNOBJ)
	ar -csru $@ $(LALNOBJ)
$(LALND)/%.o: $(LALND)/%.c
	gcc -c $(CFLAGS) -I$(LUTILS_DIR) -I$(INCLUDE)/klib $< -o $@
clean_aln:
	rm -f $(LALND)/*.o lib/aln/libaln.a

src/pileup.o: src/pileup.c
	gcc -c $(CFLAGS) -I$(LHTSLIB_INCLUDE) -I$(LUTILS_DIR) $< -o $@

src/tview.o: src/tview.c
	$(CC) -c $(CFLAGS) -I$(LHTSLIB_INCLUDE) -I$(LUTILS_DIR) $< -o $@

src/markdup.o: src/markdup.c
	$(CC) -c $(CFLAGS) -I$(LHTSLIB_INCLUDE) -I$(LUTILS_DIR) $< -o $@

src/ndr.o: src/ndr.c
	$(CC) -c $(CFLAGS) -I$(LUTILS_DIR) -I$(LKLIB_DIR) $< -o $@

src/vcf2bed.o: src/vcf2bed.c
	$(CC) -c $(CFLAGS) -I$(LUTILS_DIR) -I$(LKLIB_DIR) $< -o $@

src/epiread.o: src/epiread.c
	$(CC) -c $(CFLAGS) -I$(LHTSLIB_INCLUDE) -I$(LUTILS_DIR) $< -o $@

src/asm_pairwise.o: src/asm_pairwise.c
	$(CC) -c $(CFLAGS) -I$(LUTILS_DIR) -I$(LSGSL_DIR) $< -o $@

src/bsstrand.o: src/bsstrand.c
	$(CC) -c $(CFLAGS) -I$(LUTILS_DIR) -I$(LHTSLIB_INCLUDE) $< -o $@

src/cinread.o: src/cinread.c
	$(CC) -c $(CFLAGS) -I$(LUTILS_DIR) -I$(LHTSLIB_INCLUDE) $< -o $@

src/bsconv.o: src/bsconv.c
	$(CC) -c $(CFLAGS) -I$(LUTILS_DIR) -I$(LHTSLIB_INCLUDE) $< -o $@

src/mergecg.o: src/mergecg.c
	$(CC) -c $(CFLAGS) -I$(LUTILS_DIR) -I$(LHTSLIB_INCLUDE) $< -o $@

src/bamfilter.o: src/bamfilter.c
	$(CC) -c $(CFLAGS) -I$(LUTILS_DIR) -I$(LHTSLIB_INCLUDE) $< -o $@

# ####### general #######

# VPATH = src

# src/%.o: %.c
# 	gcc -c $(CFLAGS) $< -o $@

####### clean #######

## clean just src
.PHONY: clean
clean :
	rm -f src/*.o biscuit

## clean src and library objects
purge : clean
	make -C $(LKLIB_DIR) purge
	make -C $(LHTSLIB_DIR) clean
	make -C $(LUTILS_DIR) purge
	make -C $(LSGSL_DIR) purge
	rm -f $(LALND)/*.o $(LALND)/*.a
	rm -f biscuit

## clean to make a release zip
.PHONY: release
release:
	rm -rf release-source.zip biscuit-release
	git clone --recursive . biscuit-release
	make -C biscuit-release cleanse
	zip -r release-source.zip biscuit-release
	rm -rf biscuit-release

# removes git history, for release internal use
cleanse : purge
	rm -f **/*.o .travis.yml .gitmodules .gitignore
	rm -rf .git $(LKLIB_DIR)/.git $(LHTSLIB_DIR)/.git $(LUTILS_DIR)/.git $(LSGSL_DIR)/.git docker




