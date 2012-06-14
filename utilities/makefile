compiler = clang
compile_options = -O3 -msse2 -Wall -Wextra -D_FILE_OFFSET_BITS=64 # --short-enums

# Debugging
#compiler = gcc
#compile_options = -O0 -g -ftree-vectorize -msse2 -msse -mmmx --unroll-loops --fast-math -Wall -Wextra -D_FILE_OFFSET_BITS=64 #--short-enums 

all: verify_mapped_read_locations

verify_mapped_read_locations : verify_mapped_read_locations.c
	$(compiler) $(compile_options) -o verify_mapped_read_locations verify_mapped_read_locations.c \
	../src/rawread.o ../src/mapped_read.o ../src/fragment_length.o ../src/trace.o \
	../src/quality.o ../src/genome.o  ../src/sam.o  ../src/pseudo_location.o \
	../src/index_genome.o ../src/sequences_node.o ../src/dna_sequence.o ../src/mapped_location.o  \
	../src/iterative_mapping.o ../src/config_parsing.o ../src/diploid_map_data.o \
	../src/candidate_mapping.o ../src/error_correction.o \
	-lm -pthread
	cp verify_mapped_read_locations ../bin/

.PHONY : clean
clean :
	-rm verify_mapped_read_locations
