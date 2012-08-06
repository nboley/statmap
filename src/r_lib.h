#include <R.h>

#define R_HOME "/usr/lib/R/"

int 
load_statmap_source()
{
    char* r_code_fname = "/home/nboley/Desktop/statmap/R_lib/error_model.R";
    eval( lang2(install("source"), mkString(r_code_fname)), R_GlobalEnv );
    return 0;
}

void
init_R()
{
    fprintf( stderr, "NOTICE      :  Setting R_HOME to '%s'\n", R_HOME );
    setenv("R_HOME", R_HOME, false);

    char *argv[] = {"Rstatmap", "--gui=none", "--vanilla", "--slave"};
    int argc = sizeof(argv)/sizeof(argv[0]);
    
    Rf_initEmbeddedR(argc, argv);
    
    load_statmap_source();
    return;
}

void
end_R()
{
    Rf_endEmbeddedR(0);
    return;
}
