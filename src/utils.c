#include <R.h>
#include <Rinternals.h>

#include <stdbool.h>
#include <float.h>

#include "matdat.h"



//  sn      pointer to n

//  returns matrix, with 2 columns, of all pairs, in order matching combinations(n,2)
//  For example, allpairs(5) returns these 5*(5-1)/2 = 10 pairs:
//       [,1] [,2]
//  [1,]    1    2
//  [2,]    1    3
//  [3,]    1    4
//  [4,]    1    5
//  [5,]    2    3
//  [6,]    2    4
//  [7,]    2    5
//  [8,]    3    4
//  [9,]    3    5
// [10,]    4    5

SEXP
allpairs( SEXP sn )
    {
#if 1
    int n = *INTEGER(sn);
#else
    int n = Rf_length(sground) ;

    int *ground = INTEGER(sground);

    if( n == 1 )
        {
        //  big simplification
        n       = ground[0] ;
        ground  = NULL ;
        }
#endif

    int rowsout = ( n*(n-1) )/2 ;

    SEXP    out = PROTECT( Rf_allocMatrix(INTSXP,rowsout,2) );
    int     *matout = INTEGER(out);

    int     count   = 0;
    for( int j1=1 ; j1<n ; j1++ )
        {
        for( int j2=j1+1 ; j2<=n ; j2++, count++ )
            {
            int *p1 = matout + count ;
            int *p2 = p1 + rowsout ;

            p1[0] =  j1 ;
            p2[0] =  j2 ;
            }
        }


#if 0
    if( ground )
        {
        int m = Rf_length(out) ;

        for( int i=0 ; i<m ; i++ )
            matout[i]   = ground[ matout[i]-1 ];
        }
#endif

    UNPROTECT(1);

    return( out );
    }





//  sx      pointer to matrix of doubles,  prevalidated
//  smargin pointer to an integer, 1 or 2, prevalidated
//

//  returns a list with items
//      *)  vector of integers with col index of maximum in each row
//      *)  vector with the max values in each row

SEXP
whichMaxMatrix( SEXP sx, SEXP smargin )
    {
    matdat  md = extractmatdat( sx, smargin );
    if( md.mat == NULL )
        {
        return(R_NilValue);
        }

    SEXP    sidx    = PROTECT( Rf_allocVector(INTSXP,md.nVec) );
    int     *idxvec = INTEGER(sidx);

    SEXP    sval    = PROTECT( Rf_allocVector(REALSXP,md.nVec) );
    double  *valvec = REAL(sval);

    for( int j=0 ; j<md.nVec ; j++ )
        {
        int     idx = 0;
        double  *vec = md.mat + j*md.vecStep ;

        double  vmax = -FLT_MAX ;

        for( int i=0 ; i<md.vecLen ; i++ )
            {
            double  v = vec[i*md.eltStep] ;

            if( R_IsNA(v) )
                {
                vmax    = R_NaReal;
                idx     = R_NaInt;
                break;
                }

            if( vmax < v )
                {
                vmax    = v ;
                idx     = i + 1;      //  convert from 0-based to 1-based
                }
            }

        idxvec[j]   = idx ;
        valvec[j]   = vmax ;
        }

    SEXP    out = PROTECT( Rf_allocVector(VECSXP,2) );

    SET_VECTOR_ELT( out, 0, sidx );
    SET_VECTOR_ELT( out, 1, sval );

    UNPROTECT(3);   //  sidx and sval and variable out

    return(out);
    }



//  sx      pointer to matrix of doubles,  prevalidated
//  smargin pointer to an integer, 1 or 2, prevalidated
//
//  returns cumsum of rows or columns, depending on smargin
//
//  this newer version uses "long double" as the accumulator
//  the time is the same and accuracy is much better

SEXP
cumsumMatrix( SEXP sx, SEXP smargin )
    {
    matdat  md = extractmatdat( sx, smargin );
    if( md.mat == NULL )
        {
        return(R_NilValue);
        }

    int nrow    = md.dim[0] ;
    int ncol    = md.dim[1] ;

    SEXP    out = PROTECT( Rf_allocMatrix(REALSXP,nrow,ncol) );
    double  *matout = REAL(out);

    //  outer loop over elements
    for( int i=0 ; i<md.vecLen ; i++ )
        {
        const   double  *vecin  = md.mat + i*md.eltStep ;
        double          *vecout = matout + i*md.eltStep ;

        long double accum=0;

        //  inner loop over vectors
        for( int j=0,k=0 ; j<md.nVec ; j++, k+=md.vecStep )
            {
            accum   += (long double) vecin[k] ;

            vecout[k] = accum ;
            }
        }

    UNPROTECT(1);   //  variable out

    return(out);
    }




//  smat    pointer to matrix of doubles,  prevalidated
//  svec    a vector of doubles, with the same length as either rows or columns of smat, depending on smargin
//  smargin pointer to an integer, 1 or 2, prevalidated
//
//  adds svec to every row or column of smat in place
//  returns TRUE or NULL

SEXP
plusEqual( SEXP smat, SEXP svec, SEXP smargin )
    {
    matdat  md = extractmatdat( smat, smargin );
    if( md.mat == NULL )
        {
        return(R_NilValue);
        }

    if( Rf_length(svec) != md.vecLen )
        {
        Rprintf( "plusEqual().  %d != %d.\n", Rf_length(svec), md.vecLen );
        return(R_NilValue);
        }

    const   double  *vec = REAL(svec);

    for( int i=0 ; i<md.nVec ; i++ )
        {
        double  *v  = md.mat + i*md.vecStep ;

        for( int j=0 ; j<md.vecLen ; j++ )
            v[ j*md.eltStep ] += vec[j] ;   // this line is the crux
        }

	SEXP    out = PROTECT( Rf_allocVector(LGLSXP,1) );

    *(LOGICAL(out)) = 1 ;

    UNPROTECT(1);

    return(out);
    }



//  smat    pointer to matrix of doubles,  prevalidated
//  svec    a vector of doubles, with the same length as either rows or columns of smat, depending on smargin
//  smargin pointer to an integer, 1 or 2, prevalidated
//
//  multiplies svec times every row or column of smat in place
//  returns TRUE or NULL

SEXP
timesEqual( SEXP smat, SEXP svec, SEXP smargin )
    {
    matdat  md = extractmatdat( smat, smargin );
    if( md.mat == NULL )
        {
        return(R_NilValue);
        }

    if( Rf_length(svec) != md.vecLen )
        {
        Rprintf( "timesEqual().  %d != %d.\n", Rf_length(svec), md.vecLen );
        return(R_NilValue);
        }

    const   double  *vec = REAL(svec);

    for( int i=0 ; i<md.nVec ; i++ )
        {
        double  *v  = md.mat + i*md.vecStep ;

        for( int j=0 ; j<md.vecLen ; j++ )
            v[ j*md.eltStep ] *= vec[j] ;   // this line is the crux
        }

	SEXP    out = PROTECT( Rf_allocVector(LGLSXP,1) );

    *(LOGICAL(out)) = 1 ;

    UNPROTECT(1);

    return(out);
    }


    
    
    
    
//  smat    pointer to matrix of doubles,  prevalidated
//  svec    a vector of doubles, with the same length as either rows or columns of smat, depending on smargin
//  smargin pointer to an integer, 1 or 2, prevalidated
//
//  adds svec to every row or column of smat and returns the new matrix
//  returns TRUE or FALSE


SEXP
sumMatVec( SEXP smat, SEXP svec, SEXP smargin )
    {
    matdat  md = extractmatdat( smat, smargin );
    if( md.mat == NULL )
        {
        return(R_NilValue);
        }

    if( Rf_length(svec) != md.vecLen )
        {
        Rprintf( "sumMatVec().  %d != %d.\n", Rf_length(svec), md.vecLen );
        return(R_NilValue);
        }

    const   double  *vec = REAL(svec);

    SEXP    out = PROTECT( Rf_allocMatrix(REALSXP,md.dim[0],md.dim[1]) );
    double  *matout = REAL(out);

    for( int i=0 ; i<md.nVec ; i++ )
        {
        const   double  *vin  = md.mat + i*md.vecStep ;

        double  *vout = matout + i*md.vecStep ;

        for( int j=0 ; j<md.vecLen ; j++ )
            vout[ j*md.eltStep ]   = vin[ j*md.eltStep ] + vec[j] ;
        }

    UNPROTECT(1);

    return(out);
    }


//  smat    N x M matrix
//
//  returns rbind(mat,-mat),  but faster

SEXP
extend_antipodal( SEXP smat )
    {
    const int *dim  = INTEGER(Rf_getAttrib(smat, R_DimSymbol));

    int n   = dim[0] ;
    int m   = dim[1] ;

    const   double  *mat = REAL(smat);

    //  allocate matrix with twice the number of rows
    SEXP    out = PROTECT( Rf_allocMatrix(REALSXP,2*n,m) );
    double  *matout = REAL(out);

    for( int j=0 ; j<m ; j++ )
        {
        //  top half
        const double    *vecin  = mat + j*n ;
        double          *vecout = matout + j*2*n ;

        memcpy( vecout, vecin, n * sizeof(*vecout) );

        //  bottom half - change the sign
        vecout += n ;

        for( int i=0 ; i<n ; i++ )
            vecout[i] = -vecin[i] ;
        }

    UNPROTECT(1);   //  out

    return(out);
    }

    
////////////////            deadwood below      //////////////////////////

//  sx      pointer to matrix of doubles,  prevalidated
//  smargin pointer to an integer, 1 or 2, prevalidated
//
//  returns cumsum of rows or columns, depending on smargin

SEXP
cumsumMatrix_old( SEXP sx, SEXP smargin )
    {
    matdat  md = extractmatdat( sx, smargin );
    if( md.mat == NULL )
        {
        return(R_NilValue);
        }

    int nrow    = md.dim[0] ;
    int ncol    = md.dim[1] ;

    SEXP    out = PROTECT( Rf_allocMatrix(REALSXP,nrow,ncol) );
    double  *matout = REAL(out);

    //  the first vector is just a copy of the input
    for( int i=0, k=0 ; i<md.vecLen ; i++, k+=md.eltStep )
        matout[k] = md.mat[k] ;

    //  for( int i=0 ; i<md.vecLen ; i++ )
    //      matout[i*md.eltStep]    = md.mat[i*md.eltStep] ;

    for( int j=1 ; j<md.nVec ; j++ )
        {
        const   double  *vecin      = md.mat + j*md.vecStep ;
        double          *vecout     = matout + j*md.vecStep ;
        const   double  *vecoutprev = matout + (j-1)*md.vecStep ;

        for( int i=0, k=0 ; i<md.vecLen ; i++, k+=md.eltStep )
            {
            vecout[k] = vecoutprev[k] + vecin[k];
            }
        }

    UNPROTECT(1);   //  variable out

    return(out);
    }

