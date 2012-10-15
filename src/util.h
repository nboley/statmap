#ifndef UTIL_H
#define UTIL_H

void
safe_mkdir(char* dir);

void
safe_copy_into_output_directory( char* fname, char* output_dir, char* output_fname );

void
safe_link_into_output_directory( char* fname, char* output_dir, char* output_fname );

FILE* 
open_check_error( char* fname, char* file_mode );

void*
safe_malloc( size_t num_bytes );

void*
safe_realloc( void* ptr, size_t size );

enum bool
file_is_empty( FILE* fp );

int
cmp_ints( const int* i1,
	  const int* i2 );

#endif /* UTIL_H */
