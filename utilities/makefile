compiler = clang
#compile_options = -fast
compile_options = -O0 -g -Wall -Wextra -D_FILE_OFFSET_BITS=64 # --short-enums

all: build_index aggregate_over_wiggles mapped_reads_into_sam sample_mapping \
	estimate_read_density call_peaks \
	bootstrap_wiggles_from_mapped_reads build_random_start_from_mapped_reads \
	aggregate_over_traces convert_trace_into_wiggle

build_index : build_index.c
	$(compiler) $(compile_options) -o build_index build_index.c \
	../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/index_genome.o ../src/sequences_node.o ../src/mapped_location.o   \
	../src/genome.o ../src/pseudo_location.o ../src/snp.o \
	../src/dna_sequence.o ../src/quality.o \
	-lm -pthread
	cp build_index ../bin/

aggregate_over_wiggles : aggregate_over_wiggles.c
	$(compiler) $(compile_options) -o aggregate_over_wiggles aggregate_over_wiggles.c \
	../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/index_genome.o ../src/sequences_node.o ../src/mapped_location.o   \
	../src/genome.o ../src/pseudo_location.o ../src/snp.o \
	../src/dna_sequence.o ../src/quality.o \
	../src/wiggle.o \
	-lm -pthread
	cp aggregate_over_wiggles ../bin/

aggregate_over_traces : aggregate_over_traces.c
	$(compiler) $(compile_options) -o aggregate_over_traces aggregate_over_traces.c ../src/trace.o \
	-pthread
	cp aggregate_over_traces ../bin/

convert_trace_into_wiggle : convert_trace_into_wiggle.c
	$(compiler) $(compile_options) -o convert_trace_into_wiggle convert_trace_into_wiggle.c \
	../src/trace.o -pthread
	cp aggregate_over_traces ../bin/


mapped_reads_into_sam : mapped_reads_into_sam.c
	$(compiler) $(compile_options) -o mapped_reads_into_sam write_sam_from_mapped_reads_db.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o  ../src/pseudo_location.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	../src/iterative_mapping.o ../src/config_parsing.o \
	-lm -pthread
	cp mapped_reads_into_sam ../bin/

mapped_reads_into_sam.c :

estimate_read_density : estimate_read_density.c
	$(compiler) $(compile_options) -o estimate_read_density estimate_read_density.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o ../src/wiggle.o \
	 ../src/pseudo_location.o ../src/iterative_mapping.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	-lm -pthread
	cp estimate_read_density ../bin/

estimate_read_density.c :

sample_mapping : sample_mapping.c
	$(compiler) $(compile_options) -o sample_mapping sample_mapping.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o ../src/wiggle.o \
	 ../src/pseudo_location.o ../src/iterative_mapping.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	-lm -pthread
	cp sample_mapping ../bin/

bootstrap_wiggles_from_mapped_reads : bootstrap_wiggles_from_mapped_reads.c
	$(compiler) $(compile_options) -o bootstrap_wiggles_from_mapped_reads bootstrap_wiggles_from_mapped_reads.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o ../src/wiggle.o \
	 ../src/pseudo_location.o ../src/iterative_mapping.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	-lm -pthread
	cp bootstrap_wiggles_from_mapped_reads ../bin/

bootstrap_wiggles_from_mapped_reads.c :

build_random_start_from_mapped_reads : build_random_start_from_mapped_reads.c
	$(compiler) $(compile_options) -o build_random_start_from_mapped_reads build_random_start_from_mapped_reads.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o ../src/wiggle.o \
	 ../src/pseudo_location.o ../src/iterative_mapping.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	-lm -pthread
	cp build_random_start_from_mapped_reads ../bin/

build_random_start_from_mapped_reads.c :

call_peaks : call_peaks.c
	$(compiler) $(compile_options) -o call_peaks call_peaks.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o ../src/wiggle.o \
	 ../src/pseudo_location.o ../src/iterative_mapping.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o \
	../src/mapped_location.o  ../src/peak_calling.o \
	-lm -pthread
	cp call_peaks ../bin/

call_peaks.c :

#wiggle.o : ../src/wiggle.c
#	$(compiler) $(compile_options) -c ../src/wiggle.c ../src/trace.c ../src/mapped_read.c

.PHONY : clean
clean :
	-rm aggregate_over_wiggles
	-rm aggregate_over_traces
	-rm convert_trace_into_wiggle
	-rm mapped_reads_into_sam
	-rm estimate_read_density
	-rm build_index
	-rm call_peaks
	-rm sample_mapping
	-rm bootstrap_wiggles_from_mapped_reads
	-rm build_random_start_from_mapped_reads