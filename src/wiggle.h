
#ifndef WIGGLE_H
#define WIGGLE_H

#include "trace.h"

extern void
aggregate_over_wiggles(
    FILE** wig_fps,
    int num_wigs,
    FILE* ofp
);

#if 0
/* TODO - implement this */
extern void
load_wiggle_into_trace(
    /* the trace to init and store the data into */
    struct trace_t* trace,
    /* contains the chr names and lengths */
    struct genome_data* genome;

    /* fp's for each wiggle - each must contain exactly one track */
    FILE** ifs,
    /* the number of wiggle files */
    int num_wigs
);
#endif

extern void
write_wiggle_from_trace( 
    struct trace_t* traces,
    
    /* These are usually chr names */
    char** scaffold_names,
    /* The names of the various tracks */
    /* Use null for the indexes */
    char** track_names,
    
    const char* output_fname,                           
    const double filter_threshold 
);

#endif // #define WIGGLE_H
