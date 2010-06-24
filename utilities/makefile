compile_options = -O0 -g -Wall -Wextra -D_FILE_OFFSET_BITS=64

all: aggregate_over_wiggles mapped_reads_into_sam

aggregate_over_wiggles : aggregate_over_wiggles.c wiggle.o
	gcc $(compile_options) -o aggregate_over_wiggles aggregate_over_wiggles.c wiggle.o

mapped_reads_into_sam : 
	gcc -O0 -g -o mapped_reads_into_sam write_sam_from_mapped_reads_db.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o ../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o \
	-lm

wiggle.o : ../src/wiggle.c
	gcc $(compile_options) -c ../src/wiggle.c

.PHONY : clean
clean :
	-rm aggregate_over_wiggles wiggle.o
	-rm mapped_reads_into_sam