all: aggregate_over_wiggles

aggregate_over_wiggles : aggregate_over_wiggles.c wiggle.o
	gcc -o aggregate_over_wiggles aggregate_over_wiggles.c ../src/wiggle.o

wiggle.o : ../src/wiggle.c
	gcc -c ../src/wiggle.c

.PHONY : clean
clean :
	-rm aggregate_over_wiggles wiggle.o