compiler = clang
#compile_options = -fast
compile_options = -O0 -g -Wall -Wextra -D_FILE_OFFSET_BITS=64 # --short-enums

all: build_index aggregate_over_wiggles mapped_reads_into_sam sample_mapping \
	call_peaks aggregate_over_traces convert_trace_into_wiggle \
	build_marginal_mappings_trace_from_mpd_read_db verify_mapped_read_locations

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

verify_mapped_read_locations : verify_mapped_read_locations.c
	$(compiler) $(compile_options) -o verify_mapped_read_locations verify_mapped_read_locations.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o  ../src/pseudo_location.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	../src/iterative_mapping.o ../src/config_parsing.o \
	-lm -pthread
	cp verify_mapped_read_locations ../bin/

mapped_reads_into_sam.c :


mapped_reads_into_sam : mapped_reads_into_sam.c
	$(compiler) $(compile_options) -o mapped_reads_into_sam write_sam_from_mapped_reads_db.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o  ../src/pseudo_location.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	../src/iterative_mapping.o ../src/config_parsing.o \
	-lm -pthread
	cp mapped_reads_into_sam ../bin/

mapped_reads_into_sam.c :

build_marginal_mappings_trace_from_mpd_read_db : build_marginal_mappings_trace_from_mpd_read_db.c
	$(compiler) $(compile_options) -o build_marginal_mappings_trace_from_mpd_read_db \
	build_marginal_mappings_trace_from_mpd_read_db.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o  ../src/pseudo_location.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	../src/iterative_mapping.o ../src/config_parsing.o \
	-lm -lpthread
	cp build_marginal_mappings_trace_from_mpd_read_db ../bin/

build_marginal_mappings_trace_from_mpd_read_db.c :

sample_mapping : sample_mapping.c
	$(compiler) $(compile_options) -o sample_mapping sample_mapping.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/snp.o ../src/sam.o ../src/wiggle.o \
	 ../src/pseudo_location.o ../src/iterative_mapping.o ../src/config_parsing.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	-lm -pthread
	cp sample_mapping ../bin/

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
	-rm build_index
	-rm call_peaks
	-rm sample_mapping
	-rm build_marginal_mappings_trace_from_mpd_read_db
	-rm verify_mapped_read_locations