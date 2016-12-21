CC     = gcc
AR     = ar
CFLAGS = -g -Wall -O2

LIBUTILS_OBJS = \
	encode.o \
	stats.o \
	chisq.o \
	wzhmm.o

.c.o :
	$(CC) -c $(CFLAGS) $< -o $@

libutils.a: $(LIBUTILS_OBJS)
	@-rm -f $@
	$(AR) -csr $@ $(LIBUTILS_OBJS)

clean:
	rm -f *.o *.a
