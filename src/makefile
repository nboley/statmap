compile_options = -O3  --unroll-loops --fast-math -Wall -Wextra --short-enums -D_FILE_OFFSET_BITS=64 

objects = statmap.o genome.o index_genome.o rawread.o sequences_node.o \
	  dna_sequence.o quality.o mapped_location.o find_candidate_mappings.o \
	  iterative_mapping.o snp.o trace.o candidate_mapping.o mapped_read.o \
	  fragment_length.o wiggle.o sam.o pseudo_location.o

all: statmap

statmap : $(objects)
	gcc -o statmap \
	$(objects) \
	$(compile_options) \
	-lm -lpthread
	cp statmap ../bin/

statmap.o : statmap.c statmap.h find_candidate_mappings.h \
	index_genome.h quality.h config.h snp.h fragment_length.h sam.h
	gcc $(compile_options) -c statmap.c

snp.o : snp.c snp.h genome.h config.h
	gcc $(compile_options) -c snp.c

find_candidate_mappings.o : find_candidate_mappings.c find_candidate_mappings.h \
			    quality.h index_genome.h snp.h
	gcc $(compile_options) -c find_candidate_mappings.c

find_candidate_mappings.h : genome.h rawread.h 

#############################################################################
#
# Iterative Mapping
#

iterative_mapping.o : iterative_mapping.c iterative_mapping.h mapped_location.h statmap.h trace.h mapped_location.h snp.h genome.h fragment_length.h
	gcc $(compile_options) -c iterative_mapping.c

iterative_mapping.h : config.h genome.h trace.h
#############################################################################
#
# Data type objects
#
# 0) raw_read - the totally unmapped reads
# 1) mapped_locations - data returned byt the index search methods
# 2) candidate_mappings - mapped_locations with meta information
# 3) mapped_read - the joined cnadidate mappings
#

# utility objects

rawread.o : rawread.c config.h rawread.h quality.h dna_sequence.h
	gcc $(compile_options) -c rawread.c

rawread.h : config.h

mapped_location.o : mapped_location.c rawread.o mapped_location.h genome.h
	gcc $(compile_options) -c mapped_location.c

mapped_location.h : 

mapped_read.o : mapped_read.c mapped_read.h trace.h fragment_length.h fragment_length.o iterative_mapping.h 
	gcc $(compile_options) -c mapped_read.c

mapped_read.h : iterative_mapping.h trace.h

candidate_mapping.o : candidate_mapping.c candidate_mapping.h
	gcc $(compile_options) -c candidate_mapping.c

candidate_mapping.h : rawread.h


#############################################################################
#
# Index ( and related ) Files
#

index_genome.o : index_genome.c index_genome.h quality.h
	gcc $(compile_options) -c index_genome.c

index_genome.h : config.h genome.h sequences_node.h mapped_location.h

sequences_node.o : sequences_node.c sequences_node.h quality.h pseudo_location.h
	gcc $(compile_options) -c sequences_node.c

sequences_node.h : config.h mapped_location.h

genome.o : genome.c genome.h index_genome.h pseudo_location.h
	gcc $(compile_options) -c genome.c

genome.h : 

pseudo_location.o : pseudo_location.c pseudo_location.h
	gcc $(compile_options) -c pseudo_location.c

pseudo_location.h : config.h

#############################################################################
#
# Utility objects
#

sam.o: sam.c sam.h genome.o mapped_read.o rawread.o
	gcc $(compile_options) -c sam.c

sam.h : 

wiggle.o: wiggle.c wiggle.h
	gcc $(compile_options) -c wiggle.c

wiggle.h: trace.h

# Store traces across the gnome ( like wiggles )
trace.o: trace.c trace.h
	gcc $(compile_options) -c trace.c

trace.h:

fragment_length.o: fragment_length.c fragment_length.h
	gcc $(compile_options) -c fragment_length.c

fragment_length.h: config.h


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
