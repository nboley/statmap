compile_options = -O0 -g -Wall -Wextra -D_FILE_OFFSET_BITS=64

all: aggregate_over_wiggles

aggregate_over_wiggles : aggregate_over_wiggles.c wiggle.o
	gcc $(compile_options) -o aggregate_over_wiggles aggregate_over_wiggles.c wiggle.o

wiggle.o : ../src/wiggle.c
	gcc $(compile_options) -c ../src/wiggle.c

.PHONY : clean
clean :
	-rm aggregate_over_wiggles wiggle.o