compile_options = -O0 -g -Wall -Wextra -D_FILE_OFFSET_BITS=64

all: aggregate_over_wiggles mapped_reads_into_sam mapped_reads_into_wiggle

aggregate_over_wiggles : aggregate_over_wiggles.c
	gcc $(compile_options) -o aggregate_over_wiggles aggregate_over_wiggles.c \
	../src/wiggle.o ../src/trace.o ../src/mapped_read.o ../src/fragment_length.o ../src/pseudo_location.o \
	-lm -pthread

mapped_reads_into_sam : 
	gcc -O0 -g -o mapped_reads_into_sam write_sam_from_mapped_reads_db.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o  ../src/pseudo_location.o \
	-lm -pthread

mapped_reads_into_wiggle : 
	gcc -O0 -g -o mapped_reads_into_wiggle write_wiggle_from_mapped_reads_db.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o ../src/wiggle.o \
	 ../src/pseudo_location.o \
	-lm -pthread

#wiggle.o : ../src/wiggle.c
#	gcc $(compile_options) -c ../src/wiggle.c ../src/trace.c ../src/mapped_read.c

.PHONY : clean
clean :
	-rm aggregate_over_wiggles
	-rm mapped_reads_into_sam
	-rm mapped_reads_into_wiggle