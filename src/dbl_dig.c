#include <R.h>
#include <Rinternals.h>
#include <float.h>   /* DBL_DIG */

//  #include "config.h"         // this file is created by configure.win or configure

SEXP dbl_dig(void)
{
    SEXP out = PROTECT( Rf_allocVector(INTSXP, 1) );
    INTEGER(out)[0] = DBL_DIG;
    UNPROTECT(1);
    return out;
}
