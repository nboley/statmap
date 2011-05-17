#define max_num_qual_scores 256

struct error_data_t {
    int block_size;

    int num_unique_reads;

    int max_read_length;
    int* position_mismatch_cnts;
    
    int qual_score_cnts[max_num_qual_scores];
    int qual_score_mismatch_cnts[max_num_qual_scores];
};

void
init_error_data( 
    struct error_data_t** data,
    int max_read_length
);

void 
free_error_data( struct error_data_t* data );

void
update_error_data(
    struct error_data_t* data,
    char* genome_seq,
    char* read,
    char* error_str,
    int length
);

void
fprintf_error_data( FILE* stream, struct error_data_t* data );
