#include <math.h>

#define NDEBUG 
#include <assert.h>

#include "../src/mapped_read.h"

#define NUM_STEPS 1000

int main()
{
    int i;
    for( i = -NUM_STEPS; i < 0; i++ )
    {
        float val = pow(10, ((double)(i*70.0))/(2*NUM_STEPS) );
        // printf( "%e\t%e\n", ((double)(i*37.0))/(2*NUM_STEPS), val );
        float conv_value = ML_PRB_TYPE_to_float( ML_PRB_TYPE_from_float( val ) );
        printf("%e\t%e\t%e\t%e\n", val, conv_value, val-conv_value, (val-conv_value)/val );
    }
}
