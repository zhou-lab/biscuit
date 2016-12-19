CC     = gcc
AR     = ar
CFLAGS = -g -Wall -O2

KLIB_OBJS = \
	bgzf.o \
	kexpr.o \
	khmm.o \
	kmath.o \
	knetfile.o \
	knhx.o \
	kopen.o \
	ksa.o \
	kson.o \
	kstring.o \
	ksw.o \
	kthread.o \
	kurl.o

.c.o :
	$(CC) -c $(CFLAGS) $< -o $@

klib.a: $(KLIB_OBJS)
	@-rm -f $@
	$(AR) -csr $@ $(KLIB_OBJS)

clean:
	rm -f *.o *.a
