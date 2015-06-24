CFLAGS=-W -Wall -finline-functions -fPIC -std=gnu99

OS := $(shell uname)
ifeq ($(OS),  Darwin)
  CFLAGS += -Wno-unused-function
endif

INCLUDE = include

LSAM0119D = lib/libsamtools-0.1.19
LSAM0119 = $(LSAM0119D)/libsam.a

KLIBD = lib/klib
KLIB = lib/klib/klib.a

LUTILSD = lib/utils
LUTILS = lib/utils/libutils.a

PROG = pileup_cytosine
# PROG = hemifinder correct_bsstrand get_unmapped sample_trinuc

release : CFLAGS += -O3
release : $(PROG)

debug : CFLAGS += -g
debug : $(PROG)

$(LSAM0119) :
	make -C $(LSAM0119D) libsam.a

.PHONY: klib
klib: $(KLIB)
KLIBOBJ = $(KLIBD)/kstring.o
$(KLIB): $(KLIBOBJ)
	ar -csru $@ $(KLIBOBJ)
$(KLIBD)/%.o: $(KLIBD)/%.c
	gcc -c $(CFLAGS) -I$(INCLUDE) $< -o $@
clean_klib:
	rm -f $(KLIBD)/*.o $(KLIB)

.PHONY: utils
utils: $(LUTILS)
LUTILSOBJ = $(LUTILSD)/encode.o
$(LUTILS): $(LUTILSOBJ)
	ar -csru $@ $(LUTILSOBJ)
$(LUTILSD)/%.o: $(LUTILSD)/%.c
	gcc -c $(CFLAGS) -I$(INCLUDE) $< -o $@
clean_utils:
	rm -f $(LUTILSD)/*.o $(LUTILS)

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

.PHONY: pileup_cytosine
pileup_cytosine : bin/pileup_cytosine
bin/pileup_cytosine: $(LSAM0119) $(LUTILS) src/pileup_cytosine/pileup_cytosine.c
	gcc $(CFLAGS) -o $@ -I$(LSAM0119D) -I$(INCLUDE) src/pileup_cytosine/pileup_cytosine.c $(LSAM0119) $(LUTILS) -lpthread  -lz
clean_pileup_cytosine:
	rm -f bin/pileup_cytosine

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
