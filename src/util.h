#ifndef UTIL_H
#define UTIL_H

void
safe_mkdir(char* dir);

void
safe_copy_into_output_directory( char* fname, char* output_dir, char* output_fname );

void
safe_link_into_output_directory( char* fname, char* output_dir, char* output_fname );

#endif /* UTIL_H */
