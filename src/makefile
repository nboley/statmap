CC=clang
CFLAGS=-O3 -g -msse2 -Wall -Wextra -D_FILE_OFFSET_BITS=64
objects := $(patsubst %.c,%.o,$(wildcard *.c))

# Debugging
#CC=gcc
#CFLAGS=-O0 -g -ftree-vectorize -msse2 -msse -mmmx --unroll-loops --fast-math -Wall -Wextra -D_FILE_OFFSET_BITS=64

all: statmap statmap.so

statmap: $(objects)
	$(CC) -o statmap \
	$(objects) \
	$(CFLAGS) \
	-lm -lpthread
	cp statmap ../bin/

statmap.so : $(objects)
	$(CC) $(CFLAGS) \
	-fpic -shared -Wl,-soname,libstatmap.so.1 -o libstatmap.so \
	$(wildcard *.c)

.PHONY : clean
clean :
	-rm $(objects)
	-rm statmap
	-rm libstatmap.so
