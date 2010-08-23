compile_options = -O0 -g -Wall -Wextra -D_FILE_OFFSET_BITS=64

all: build_index aggregate_over_wiggles mapped_reads_into_sam mapped_reads_into_wiggle estimate_read_density call_peaks

build_index : build_index.c
	gcc $(compile_options) -o build_index build_index.c \
	../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/index_genome.o ../src/sequences_node.o ../src/mapped_location.o   \
	../src/genome.o ../src/pseudo_location.o ../src/snp.o \
	../src/dna_sequence.o ../src/quality.o \
	-lm -pthread
	cp build_index ../bin/

aggregate_over_wiggles : aggregate_over_wiggles.c
	gcc $(compile_options) -o aggregate_over_wiggles aggregate_over_wiggles.c \
	../src/wiggle.o ../src/trace.o ../src/mapped_read.o ../src/fragment_length.o ../src/pseudo_location.o \
	-lm -pthread
	cp aggregate_over_wiggles ../bin/

mapped_reads_into_sam : 
	gcc $(compile_options) -o mapped_reads_into_sam write_sam_from_mapped_reads_db.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o  ../src/pseudo_location.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	-lm -pthread
	cp mapped_reads_into_sam ../bin/

mapped_reads_into_wiggle : 
	gcc $(compile_options) -o mapped_reads_into_wiggle write_wiggle_from_mapped_reads_db.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o ../src/wiggle.o \
	 ../src/pseudo_location.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	-lm -pthread
	cp mapped_reads_into_wiggle ../bin/

estimate_read_density : 
	gcc $(compile_options) -o estimate_read_density estimate_read_density.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o ../src/wiggle.o \
	 ../src/pseudo_location.o ../src/iterative_mapping.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	-lm -pthread
	cp estimate_read_density ../bin/

call_peaks : 
	gcc $(compile_options) -o call_peaks call_peaks.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o ../src/wiggle.o \
	 ../src/pseudo_location.o ../src/iterative_mapping.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	-lm -pthread
	cp call_peaks ../bin/


#wiggle.o : ../src/wiggle.c
#	gcc $(compile_options) -c ../src/wiggle.c ../src/trace.c ../src/mapped_read.c

.PHONY : clean
clean :
	-rm aggregate_over_wiggles
	-rm mapped_reads_into_sam
	-rm mapped_reads_into_wiggle
	-rm estimate_read_density
	-rm build_index
	-rm call_peaks