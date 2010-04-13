
compile_options = -O0 -g --unroll-loops --fast-math -Wall -Wextra --short-enums -D_FILE_OFFSET_BITS=64 

objects = statmap.o db_interface.o genome.o index_genome.o rawread.o sequences_node.o \
	  dna_sequence.o quality.o mapped_location.o find_candidate_mappings.o \
	  iterative_mapping.o snp.o

all: statmap

statmap : $(objects)
	gcc -o statmap \
	$(objects) \
	$(compile_options) \
	-lm -lpthread --openmp

statmap.o : statmap.c statmap.h find_candidate_mappings.h \
	index_genome.h db_interface.h quality.h config.h snp.h
	gcc $(compile_options) -c statmap.c

snp.o : snp.c snp.h genome.h config.h
	gcc $(compile_options) -c snp.c

find_candidate_mappings.o : find_candidate_mappings.c find_candidate_mappings.h \
			    quality.h index_genome.h snp.h
	gcc $(compile_options) -c find_candidate_mappings.c

find_candidate_mappings.h : genome.h rawread.h db_interface.h

iterative_mapping.o : iterative_mapping.c iterative_mapping.h mapped_location.h
	gcc $(compile_options) --openmp -c iterative_mapping.c

iterative_mapping.h : config.h genome.h

# utility objects
db_interface.o : db_interface.c db_interface.h mapped_location.h
	gcc $(compile_options) -c db_interface.c

db_interface.h : mapped_location.h

mapped_location.o : mapped_location.c rawread.o mapped_location.h genome.h
	gcc $(compile_options) -c mapped_location.c

mapped_location.h : config.h genome.h rawread.h

index_genome.o : index_genome.c index_genome.h quality.h
	gcc $(compile_options) -c index_genome.c

index_genome.h : config.h genome.h sequences_node.h mapped_location.h

sequences_node.o : sequences_node.c sequences_node.h quality.h
	gcc $(compile_options) -c sequences_node.c

sequences_node.h : config.h mapped_location.h

genome.o : genome.c genome.h index_genome.h
	gcc $(compile_options) -c genome.c

rawread.o : rawread.c config.h rawread.h quality.h dna_sequence.h
	gcc $(compile_options) -c rawread.c

rawread.h : config.h

quality.o : config.h quality.h
	gcc $(compile_options) -c quality.c

quality.h : dna_sequence.h

dna_sequence.o : dna_sequence.c dna_sequence.h
	gcc $(compile_options) -c dna_sequence.c

dna_sequence.h : config.h

### Clean Rules
.PHONY : clean
clean :
	-rm $(objects)
