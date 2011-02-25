/* Copyright (c) 2009-2010, Nathan Boley */

#include "config_parsing.h"

#include <stdio.h>
#include <stdlib.h>

void
write_config_file_to_disk( struct args_t* args  )
{
    /* write the arguments to a file */
    FILE* arg_fp = fopen( "config.dat", "w" );
    fwrite( &args, sizeof(struct args_t), 1, arg_fp );
    fclose( arg_fp  );
    
    return;
}

void
read_config_file_from_disk( struct args_t** args  )
{
    *args = malloc( sizeof(struct args_t)  );
    
    /* write the arguments to a file */
    FILE* arg_fp = fopen( "config.dat", "r" );
    fread( *args, sizeof(struct args_t), 1, arg_fp );
    fclose( arg_fp  );
    
    return;
}
