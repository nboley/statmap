/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <libgen.h> // for basename()

#include "../src/index_genome.h"
#include "../src/diploid_map_data.h"

int min_num_hq_bps = -1;
int num_threads = 1;

/* 
   usually, we can just let the OS cleanup the in memory index. 
   However, if we are searching for memory leaks, it may be useful 
   to manually free the index. 
*/
#define DONT_EXPLICITLY_FREE_INDEX

void usage()
{
    fprintf( stderr, "Usage: ./build_index indexed_seq_len output_filename genome.fa(s) [diploid.map(s)]\n" );
    return;
}

struct file_group {
    char* prefix;
    int num_files;
    char** filenames;
};

struct file_group_list {
    int num_groups;
    struct file_group* groups;
};

void
init_file_group (
    struct file_group* fg,
    char* prefix
)
{
    /* NOTE: assumes memory for file_group struct fg has been allocated
     * (by realloc in add_group_to_file_group_list) */

    // allocate memory for, and copy, prefix
    fg->prefix = calloc( strlen(prefix) + 1, sizeof(char) );
    strcpy( fg->prefix, prefix );

    fg->num_files = 0;
    fg->filenames = NULL;
}

void
free_file_group(
    struct file_group* fg
)
{
    // free all allocated filenames
    int i;
    for( i=0; i < fg->num_files; i++ )
    {
        free( fg->filenames[i] );
    }
    
    free( fg->filenames );
    free( fg->prefix );
}

void
init_file_group_list(
    struct file_group_list** fgl
)
{
    (*fgl) = malloc( sizeof( struct file_group_list ) );
    
    (*fgl)->num_groups = 0;
    (*fgl)->groups = NULL;
}

void
free_file_group_list(
    struct file_group_list* fgl
)
{
    // free memory allocated in file groups
    int i;
    for( i=0; i < fgl->num_groups; i++ )
    {
        free_file_group( &(fgl->groups[i]) );
    }

    // free contiguous block of file_group structs
    free( fgl->groups );
    free( fgl );
}

void
add_filename_to_file_group(
    struct file_group* fg,
    char* filename
)
{
    fg->num_files += 1;

    // Allocate a pointer to char for filename
    fg->filenames = realloc(
        fg->filenames,
        fg->num_files * sizeof(char *)
    );
    // Allocate space for the filename
    fg->filenames[fg->num_files - 1] = calloc( strlen(filename) + 1, sizeof(char) );
    // Copy filename
    strcpy( fg->filenames[fg->num_files - 1], filename );
}

void
add_group_to_file_group_list(
    char* prefix,
    struct file_group_list* fgl
)
{
    fgl->num_groups += 1;

    // Allocate memory for new file group
    fgl->groups = realloc(
        fgl->groups,
        fgl->num_groups * sizeof( struct file_group )
    );

    // Initialize new file group
    init_file_group( &(fgl->groups[fgl->num_groups - 1]), prefix );
}

void
print_file_group(
    struct file_group* fg
)
{
    // print the prefix
    printf("prefix : %s\n", fg->prefix);
    // print the filenames
    int i;
    for( i=0; i < fg->num_files; i++ )
    {
        printf("%s\n", fg->filenames[i] );
    }
}

void
print_file_group_list(
    struct file_group_list* fgl
)
{
    // print the file groups
    int i;
    for( i=0; i < fgl->num_groups; i++ )
    {
        print_file_group( &(fgl->groups[i]) );
    }
}

int
find_prefix_index_in_file_group_list(
    struct file_group_list* fgl,
    char* prefix
)
{
    int prefix_index;
    /* loop through files in group */
    for( prefix_index=0; prefix_index < fgl->num_groups; prefix_index++ )
    {
        if( 0 == strcmp( fgl->groups[prefix_index].prefix, prefix ) )
            return prefix_index;
    }

    /* if prefix not found in list of file groups, return -1 */
    return -1;
}

void
group_files_for_haploid_mapping(
    char** filenames,
    int num_files,
    struct file_group_list* fgl
)
{
    int i;
    /* loop through filenames */
    for( i = 0; i < num_files; i++ )
    {
        // default: add each file to its own group, using the full filename as the prefix
        add_group_to_file_group_list( filenames[i], fgl );
        add_filename_to_file_group( &(fgl->groups[i]), filenames[i] );
    }
}

void
group_files_for_diploid_mapping(
    char** filenames,
    int num_files,
    struct file_group_list* fgl
)
{
    int i;
    /* loop through filenames */
    for( i = 0; i < num_files; i++ )
    {
        /* find the prefix */
        /* make a copy of the original string, since basename and strtok are destructive */
        char* fncopy = calloc( strlen(filenames[i]) + 1, sizeof(char) );
        strcpy( fncopy, filenames[i] );
        char* bname = basename( fncopy );
        char* prefix = strtok( bname, "_" );

        int prefix_index = find_prefix_index_in_file_group_list( fgl, prefix );

        /* if we didn't find the prefix in the list, 
         * add a new file group to the list */
        if( prefix_index == -1 )
        {
            add_group_to_file_group_list( prefix, fgl );
            prefix_index = fgl->num_groups - 1; // set prefix_index to newly added file group
        }

        /* add the filename to the file group */
        add_filename_to_file_group( &(fgl->groups[prefix_index]), filenames[i] );

        free( fncopy );
    }
}

/* checks list of filenames for .map files
 * If it finds one, map the input as a diploid genome */
enum bool
is_diploid_genome(
    char** filenames,
    int num_files
)
{
    int i;
    for( i=0; i < num_files; i++ )
    {
        if( 0 == strcmp( ".map", filenames[i] + strlen(filenames[i]) - 4 ) )
        {
            return true;
        }
    }
    return false;
}

void
group_files(
    char** filenames,
    int num_files,
    struct file_group_list* fgl
)
{
    /* check input files to see if we're building a haploid or a diploid genome */
    enum bool is_diploid = is_diploid_genome( filenames, num_files );
    if( is_diploid )
        group_files_for_diploid_mapping( filenames, num_files, fgl );
    else
        group_files_for_haploid_mapping( filenames, num_files, fgl );
}

/* Sanity checks on input file groups */
/* Exits with fatal error if bad input is found */
void
verify_file_groups(
    struct file_group_list* fgl
)
{
    /* Iterate over groups in list */
    int i;
    for( i=0; i < fgl->num_groups; i++ )
    {
        // check number of files
        int num_files = fgl->groups[i].num_files;
        if( !(num_files == 1 || num_files == 3) )
        {
            fprintf( stderr,
                    "FATAL : File group contains %i files\n"
                    "There should be 1 (for a haploid geome), or 3 (for a diploid genome).\n",
                    num_files
            );
            exit( -1 );
        }

        // make sure diploid file groups make sense
        if( num_files == 3 )
        {
            // count types of file by suffix
            // there should be 2 .fa files and 1 .map file
            int num_fa = 0;
            int num_map = 0;
            int j;
            for( j=0; j < num_files; j++ )
            {
                char* filename = fgl->groups[i].filenames[j];
                if( 0 == strcmp( ".map", filename + strlen(filename ) - 4 ))
                    num_map++;
                else
                {
                    // Assume it's fasta
                    num_fa++;
                }
            }
            // check counts of fasta and map files
            if( !( num_fa == 2 && num_map == 1 ) )
            {
                fprintf( stderr,
                    "FATAL : Encountered a group of 3 files with %i fastas "
                    "and %i map files.\n",
                    num_fa, num_map
                );
                exit( -1 );
            }
        }
    }
}

enum CHR_SOURCE
get_chr_source_from_fasta_file(
    char* filename
)
{
    enum CHR_SOURCE chr_source;

    FILE* fp = fopen( filename, "r" );
    char id[500];
    fscanf( fp, ">%s", id );
    fclose( fp );

    // parse on underscore
    char* past_underscore = strchr( id, '_' ) + 1;

    if( past_underscore == NULL )
    {
        // fail closed
        fprintf( stderr, "FATAL : Could not determine chromosome source from .fa file %s\n", filename );
        exit( -1 );
    }

    // TODO: case insensitive? 
    if( 0 == strcmp( "paternal", past_underscore ) )
        chr_source = PATERNAL;
    else if( 0 == strcmp( "maternal", past_underscore ) )
        chr_source = MATERNAL;
    else
    {
        // fail closed
        fprintf( stderr, "FATAL : Could not determine chromosome source from .fa file %s. Found %s\n",
            filename,
            id
        );
        exit( -1 );
    }

    return chr_source;
}

void
add_file_group_to_genome(
    struct file_group* fg,
    struct genome_data* genome
)
{
    /* Is this a diploid or haploid group of genome files? */
    enum bool diploid_group;
    if( fg->num_files == 1 ) // one file (.fa) in the group => haploid
        diploid_group = false; 
    else if( fg->num_files == 3 ) // two .fa and one .map => diploid
        diploid_group = true;

    /* Loop through files in group */
    int i;
    for( i=0; i < fg->num_files; i++ )
    {
        char* filename = fg->filenames[i];
        // if file is not .map, we assume it's fasta and add it
        if( 0 != strcmp( ".map", filename + (strlen(filename) - 4)) )
        {
            /* Determine chromosome source */
            enum CHR_SOURCE chr_source; 
            if( diploid_group )
                chr_source = get_chr_source_from_fasta_file( fg->filenames[i] );
            else
                chr_source = REFERENCE;

            FILE* genome_fp = fopen( filename, "r" );
            if( NULL == genome_fp )
            {
                fprintf( stderr, "FATAL : Error opening input file %s. Source : %i", filename, chr_source );
                exit( -1 );
            }

            /* Add to genome */
            //fprintf( stdout, "NOTICE : Adding %s to genome\n", filename );
            add_chrs_from_fasta_file( genome, genome_fp, chr_source );

            //TODO: add_chrs closes input FILE*. This behavior should be consistent.
            //fclose( genome_fp );
        }
    }
}

void add_file_group_list_to_genome(
    struct file_group_list* fgl,
    struct genome_data* genome
)
{
    /* Iterate over groups in list */
    int i;
    for( i=0; i < fgl->num_groups; i++ )
    {
        add_file_group_to_genome( &(fgl->groups[i]), genome );
    }
}

int
count_map_files(
    struct file_group_list* fgl
)
{
    int num_map_files = 0;
    int i, j;
    /* loop through file group list */
    for( i=0; i < fgl->num_groups; i++ )
    {
        struct file_group* fg = &(fgl->groups[i]);
        /* loop through filenames in list */
        for( j=0; j < fg->num_files; j++ )
        {
            /* if the file is a .map file, add it to map_fnames */
            if( 0 == strcmp( ".map", fg->filenames[j] + strlen(fg->filenames[j]) - 4 ))
                num_map_files++;
        }
    }

    return num_map_files;
}

void
get_chr_name_from_fasta(
    char* filename,
    char* chr_name
)
{
    FILE* fp = fopen( filename, "r" );
    char id[500];
    fscanf( fp, ">%s", id );
    fclose( fp );

    chr_name = malloc( (strlen(id) + 1) * sizeof( char ) );
    strcpy( chr_name, id );
    /* Caller must free chr_name */
}

void
parse_diploid_file_group(
    struct file_group* fg,
    struct diploid_map_data_t* map_data,
    struct genome_data* genome
)
{
    char* paternal_fname = NULL;
    char* maternal_fname = NULL;
    char* map_fname = NULL;

    /* loop through files in group to get specific filenames */
    int i;
    for( i=0; i < fg->num_files; i++ )
    {
        char* filename = fg->filenames[i];
        if( 0 == strcmp( ".map", filename + strlen(filename) - 4 ) )
        {
            map_fname = filename;
        }
        else
        {
            enum CHR_SOURCE chr_source =  get_chr_source_from_fasta_file( filename );
            if( chr_source == PATERNAL )
                paternal_fname = filename;
            else if( chr_source == MATERNAL )
                maternal_fname = filename;
            else
            {
                fprintf( stderr, "Encountered chr of unknown origin while parsing diploid genome.\n" );
                exit(1);
            }
        }
    }

    /* get indexes of paternal and maternal chrs based on their names */
    char paternal_chr_name[500];
    FILE* pfp = fopen( paternal_fname, "r" );
    fscanf( pfp, ">%s", paternal_chr_name);
    fclose( pfp );

    char maternal_chr_name[500];
    FILE* mfp = fopen( maternal_fname, "r" );
    fscanf( mfp, ">%s", maternal_chr_name);
    fclose( mfp );

    int paternal_chr_index = find_chr_index( genome, paternal_chr_name );
    int maternal_chr_index = find_chr_index( genome, maternal_chr_name );

    /* parse map file */
    parse_map_file(
        map_fname,
        map_data,
        genome,
        paternal_chr_index,
        maternal_chr_index
    );
}

int
build_diploid_map_data(
    struct diploid_map_data_t** map_data,
    struct file_group_list* fgl,
    struct genome_data* genome
)
{
    /* count .map files in file group list */
    int num_map_files = count_map_files( fgl );
    /* If there are no .map files, return 0 for num_diploid_chrs */
    if( num_map_files == 0 )
        return 0;

    /* allocate contiguous memory; one diploid_map_data_t struct for each .map file */
    *map_data = malloc( num_map_files * sizeof( struct diploid_map_data_t ) );

    int map_data_index = 0;
    int i;
    /* loop through file group list */
    for( i=0; i < fgl->num_groups; i++ )
    {
        struct file_group* fg = &(fgl->groups[i]);
        if( is_diploid_genome( fg->filenames, fg->num_files ) )
        {
            /* parse_map_file for this diploid_map_data_t */
            parse_diploid_file_group( fg, map_data[map_data_index], genome );
            index_diploid_map_data( map_data[map_data_index] );
            map_data_index++; // increment # of map files processed
        }
    }

    return num_map_files;
}

int 
main( int argc, char** argv )
{
    if( argc < 4 )
    {
        usage();
        exit( 1 );
    }

    int indexed_seq_len = atoi( argv[1] );

    char* output_fname = argv[2];
    char index_fname[500];
    sprintf( index_fname, "%s.index", output_fname );

    /* sort input files into groups by their prefixes */
    struct file_group_list* fgl = NULL;
    init_file_group_list( &fgl );
    group_files( argv+3, argc-3, fgl );
    verify_file_groups( fgl );

    /* Initialize the genome data structure */
    struct genome_data* genome;
    init_genome( &genome );

    /*** Load the genome ***/
    add_file_group_list_to_genome( fgl, genome );

    /* Initialize the genome data structure */
    init_index( &(genome->index), indexed_seq_len );

    /* set num_diploid_chrs */
    genome->index->num_diploid_chrs = build_diploid_map_data(
            &(genome->index->map_data), fgl, genome
        );

    // DEBUG
    FILE* fp = fopen( "chr1_debug_dmap", "w" );
    write_diploid_map_data_to_file(
            genome->index->map_data,
            genome->index->num_diploid_chrs,
            fp
        );
    fclose( fp );

    /* free file group list */
    free_file_group_list( fgl );

    /* index the genome */
    index_genome( genome, indexed_seq_len );

    /* write the genome to file */
    write_genome_to_disk( genome, output_fname  );
    
    /* write the index to disk */
    build_ondisk_index( genome->index, index_fname  );

    /* Test the code - make sure the two copies are identical
    struct genome_data* genome_copy;
    load_genome_from_disk( &genome_copy, "ODGenome.bin" );
    // write the genome to file
    write_genome_to_disk( genome_copy, "ODGenome_2.bin"  );
    */

#ifdef EXPLICITLY_FREE_INDEX
    free_genome( genome );
#endif

    return 0;
}
