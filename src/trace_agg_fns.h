
#ifndef TRACE_AGG_FNS
#define TRACE_AGG_FNS

#include "../src/trace.h"

static inline TRACE_TYPE
sum( const TRACE_TYPE a, const TRACE_TYPE b )
{
    return a+b;
}

static inline TRACE_TYPE
min( const TRACE_TYPE a, const TRACE_TYPE b )
{
    if ( a < b )
        return a;
    
    return b;
}

static inline TRACE_TYPE
max( const TRACE_TYPE a, const TRACE_TYPE b )
{
    if ( a > b )
        return a;
    
    return b;
}

#endif // ifndef TRACE_AGG_FNS
