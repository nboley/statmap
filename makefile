CC=clang
CFLAGS=-O0 -g -msse2 -Wall -Wextra -D_FILE_OFFSET_BITS=64 -Wunreachable-code -Wunused
src_objects := $(patsubst %.c,%.o,$(wildcard src/*.c))

#### Global targets
all: statmap statmap.so r_error_model # verify_mapped_read_locations

clean: 
	-rm R_lib/error_model
	-rm $(src_objects)
	-rm src/statmap
	-rm src/statmap.so
	-rm utilities/verify_mapped_read_locations

check: 
	cd tests; rm -rf smo_* *.fa *.map; python tests.py && python simulate_chipseq.py;

tags:
	cd src; ctags *.c;
	cd utilities; ctags *.c *.py;
	cd python_lib; ctags *.py;
	ctags --file-scope=no -R;

### the 'R_lib' subdirectory
r_error_model: R_lib/error_model.c
	export R_HOME=/usr/lib/R
	gcc -std=gnu99 -I/usr/share/R/include -O0 -g -pipe -o R_lib/error_model R_lib/error_model.c -L/usr/lib/R/lib -lR

### the 'src' subdirectory
statmap: $(src_objects)
	$(CC) -o src/statmap \
	$(src_objects) \
	$(CFLAGS) \
	-lm -lpthread
	cp src/statmap ./bin/

statmap.so : $(src_objects)
	$(CC) $(CFLAGS) \
	-fpic -shared -Wl,-soname,libstatmap.so.1 -o src/libstatmap.so \
	$(wildcard src/*.c)

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


