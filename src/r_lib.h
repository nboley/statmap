#include <R.h>
#include "log.h"

#define R_HOME "/usr/lib/R/"

void
load_statmap_source( char* statmap_dir )
{
    char r_code_fname[512];
    strcpy( r_code_fname, statmap_dir );
    strcat( r_code_fname, "/../R_lib/error_model.R" );

    int error_occurred;
    R_tryEval( lang2(install("source"),
               mkString(r_code_fname)),
               R_GlobalEnv,
               &error_occurred );
    if( error_occurred ) {
        statmap_log( LOG_FATAL,
            "An error occurred while the R interpreter was running." );
    }
    
    return;
}

void
init_R( )
{
    statmap_log( LOG_NOTICE, "Setting R_HOME to '%s'", R_HOME );
    setenv("R_HOME", R_HOME, false);

    char *argv[] = {"Rstatmap", "--gui=none", "--vanilla", "--slave"};
    int argc = sizeof(argv)/sizeof(argv[0]);
    
    Rf_initEmbeddedR(argc, argv);
    
    return;
}

void
end_R()
{
    // Rf_endEmbeddedR(0);
    return;
}

/* get the list element named str, or return NULL */
     
SEXP getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
     
    for (R_len_t i = 0; i < length(list); i++)
        if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
            elmt = VECTOR_ELT(list, i);
            break;
        }
    return elmt;
}
