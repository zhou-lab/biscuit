CC     = gcc
AR     = ar
CFLAGS = -g -Wall -O2 -std=gnu99
# -std=c99 is important, otherwise duplicate definition for inline

SOURCES := $(wildcard **/*.c)
OBJECTS := $(patsubst %.c, %.o, $(SOURCES))

%.o : %.c
	$(CC) -c $(CFLAGS) -I. $< -o $@

# main: $(OBJECTS)
# 	@echo $(OBJECTS)

libgsl.a: $(OBJECTS)
	@-rm -f $@
	$(AR) -csr $@ $^

clean:
	rm -f $(OBJECTS)
