CC=clang
CFLAGS=-O0 -g -msse2 -Wall -Wextra -D_FILE_OFFSET_BITS=64 \
	-Wunreachable-code -Wunused \
	-I/usr/share/R/include
RLIB=/usr/lib/R/lib/libR.so
src_objects := $(patsubst %.c,%.o,$(wildcard src/*.c))

#### Global targets
all: statmap.so statmap # verify_mapped_read_locations

clean: 
	-rm $(src_objects)
	-rm src/statmap
	-rm src/statmap.so
	-rm utilities/verify_mapped_read_locations

check: 
	cd tests; rm -rf smo_* *.fa *.map; python tests.py #&& python simulate_chipseq.py;

tags:
	cd src; ctags *.c;
	cd utilities; ctags *.c *.py;
	cd python_lib; ctags *.py;

### the 'src' subdirectory
statmap: $(src_objects)
	$(CC) -o src/statmap \
	$(src_objects) $(RLIB) \
	$(CFLAGS) \
	-lm -lpthread -L/usr/lib/R/lib -lR 
	cp src/statmap ./bin/

statmap.so : $(src_objects)
	$(CC) $(CFLAGS) \
	-fpic -shared -Wl,-soname,libstatmap.so.1 -o src/libstatmap.so \
	-L/usr/lib/R/lib -lR \
	$(wildcard src/*.c) $(RLIB)

### the 'utilities' subdirectory
verify_mapped_read_locations : utilities/verify_mapped_read_locations.c
	$(compiler) $(compile_options) -o ./utilities/verify_mapped_read_locations \
	./utilities/verify_mapped_read_locations.c \
	./src/rawread.o ./src/mapped_read.o ./src/fragment_length.o ./src/trace.o \
	./src/quality.o ./src/genome.o  ./src/sam.o  ./src/pseudo_location.o \
	./src/index_genome.o ./src/sequences_node.o ./src/dna_sequence.o ./src/mapped_location.o  \
	./src/iterative_mapping.o ./src/config_parsing.o ./src/diploid_map_data.o \
	./src/candidate_mapping.o ./src/error_correction.o ./src/util.o \
	-lm -pthread
	cp utilities/verify_mapped_read_locations ./bin/

