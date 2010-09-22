
#ifndef WIGGLE_H
#define WIGGLE_H

#include "trace.h"

struct wig_line_info {
    FILE* fp;
    int file_index;
    int trace_index;
    int chr_index;
    unsigned int position;
    float value;
};

float 
wig_lines_min( const struct wig_line_info* lines, const int ub, const int num_wigs  );

float 
wig_lines_max( const struct wig_line_info* lines, const int ub, const int num_wigs  );

float 
wig_lines_sum( const struct wig_line_info* lines, const int ub, const int num_wigs  );

extern void
aggregate_over_wiggles(
    FILE** wig_fps,
    int num_wigs,
    FILE* ofp,
    float threshold,
    float agg_fn( const struct wig_line_info*, const int, const int  )
);

extern void
call_peaks_from_wiggles(
    FILE* IP_wig_fp,
    FILE* NC_wig_fp,

    struct genome_data* genome,
    struct trace_t* trace
    );


#endif // #define WIGGLE_H
