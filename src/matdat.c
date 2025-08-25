
//  #include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>

#include "matdat.h"

matdat
extractmatdat( SEXP sx, SEXP smargin )
    {
    matdat  out;

    memset( &out, 0, sizeof(out) );

    int *dim    = INTEGER(Rf_getAttrib(sx, R_DimSymbol));
    int nrow    = out.dim[0] = dim[0] ;
    int ncol    = out.dim[1] = dim[1] ;

    int type = TYPEOF(sx) ;

    if( type == REALSXP )
        out.mat     = REAL(sx) ;    //  pointer to the matrix of doubles
    else if( type == INTSXP )
        out.imat    = INTEGER(sx) ;
    else
        return( out );


    int margin  = *(INTEGER(smargin)) ;

    switch( margin )
        {
        case 1 :
            out.vecStep = 1;
            out.eltStep = nrow ;
            out.vecLen  = ncol ;
            out.nVec    = nrow ;
            break ;

        case 2 :
            out.eltStep = 1;
            out.vecStep = nrow ;
            out.vecLen  = nrow ;
            out.nVec    = ncol ;
            break ;

        default :
            out.mat = NULL;
            break;
        }

    return(out);
    }