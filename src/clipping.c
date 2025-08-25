
#include <stdbool.h>

#include "clipping.h"


//  vertex  the polygon in 3D, the last index is from 0 to n-1
//  n-1     the # of vertices
//  afun    the affine functional value at each vertex,
//          defining the positive halfspace to which to clip_poly
//  vertout the output vertex, must be allocated to hold n+1 vertices
//  nout    the # of output vertices, n+1 is max possible

bool
clip_poly( const double *vertex[3], int n, const double afun[],
            double *vertout[3], int *nout )
    {
    *nout = 0 ;

    for( int i=0 ; i<n ; i++ )
        {
        int     i_next = (i+1) % n ;

        if( afun[i] * afun[i_next] < 0 )
            {
            //  the two afun's have opposite sign
            //  the edge intersects the hyperplane
            if( 0 < afun[i] )
                {
                //  output Vi
                vertout[0][ *nout ] = vertex[0][i];
                vertout[1][ *nout ] = vertex[1][i];
                vertout[2][ *nout ] = vertex[2][i];
                (*nout)++ ;
                }

            //  compute and output the intersection
#if 1
            //  slightly modified formula for the intersection
            //  3 divisions, but coordinate which should be exactly 0, IS !
            double  denom = afun[i_next] - afun[i] ;    //  denom cannot be 0
            vertout[0][ *nout ] = (afun[i_next]*vertex[0][i] - afun[i]*vertex[0][i_next]) / denom ;
            vertout[1][ *nout ] = (afun[i_next]*vertex[1][i] - afun[i]*vertex[1][i_next]) / denom ;
            vertout[2][ *nout ] = (afun[i_next]*vertex[2][i] - afun[i]*vertex[2][i_next]) / denom ;
#else
            //  the classical formula, with only 1 division
            double  t = afun[i_next] / (afun[i_next] - afun[i]) ;
            vertout[0][ *nout ] = t*vertex[0][i]  +  (1-t)*vertex[0][i_next] ;
            vertout[1][ *nout ] = t*vertex[1][i]  +  (1-t)*vertex[1][i_next] ;
            vertout[2][ *nout ] = t*vertex[2][i]  +  (1-t)*vertex[2][i_next] ;
#endif
            (*nout)++ ;
            }
        else if( 0 <= afun[i] )
            {
            //  output Vi
            vertout[0][ *nout ] = vertex[0][i];
            vertout[1][ *nout ] = vertex[1][i];
            vertout[2][ *nout ] = vertex[2][i];
            (*nout)++ ;
            }
        }

    return( true );
    }


//  in this we only clip a quadrangle against the positive octant
//  vertout[] must be big enough to hold 7 vertices
//  if nout becomes <= 2 we have a degenerate poly and quit prematurely

bool
clipquad3D( double quad[3][4], double *vertout[3], int *nout )
    {
    //  since there are 4 vertices input, and 3 planes,
    //  there are at most 7 vertices output

    const   double  *q[3] = { quad[0], quad[1], quad[2] } ;

    double  vertout1[3][VERTSOUTMAX], vertout2[3][VERTSOUTMAX] ;

    double  *v1[3]  = { vertout1[0], vertout1[1], vertout1[2] } ;
    double  *v2[3]  = { vertout2[0], vertout2[1], vertout2[2] } ;

    bool    ok;

    *nout   = 0;

    //  plane x=0, from quad[] to vertout1[]
    ok  = clip_poly( q, 4, q[0], v1, nout );
    if( ! ok )  return(false);
    if( *nout == 0 ) return(true);

    //  plane y=0,  from vertout1[] to vertout2[]
    ok  = clip_poly( (const double **) v1, *nout, vertout1[1], v2, nout );
    if( ! ok )  return(false);
    if( *nout == 0 ) return(true);

    //  plane z=0,  from vertout2[] to final vertout[]
    ok  = clip_poly( (const double **) v2, *nout, vertout2[2], vertout, nout );

    return(ok);
    }



/*************************   testing area   *********************************/

#include <R.h>
#include <Rinternals.h>


//  these functions are only used for testing

//  smatquad        4x3 matrix

SEXP
clipquad( SEXP smatquad )
    {
    const   int *dim ;

    dim = INTEGER( Rf_getAttrib(smatquad, R_DimSymbol) );
    if( dim[0] != 4  ||  dim[1] != 3 )
        {
        Rprintf( "bad smatquad %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }

    const   double  *matquad = REAL(smatquad);

    //  const   double  *quad[3] = { matquad, matquad + 4, matquad + 2*4 };

    //  load quad
    double  quad[3][4] ;
    for( int i=0 ; i<4 ; i++ )
        for( int j=0 ; j<3 ; j++ )
            quad[j][i]  = matquad[ i + 4*j ] ;

    double  vertout[3][VERTSOUTMAX] ;

    double  *v[3]  = { vertout[0], vertout[1], vertout[2] } ;

    int nout;

    bool    ok;

    ok  = clipquad3D( quad, v, &nout ) ;
    if( ! ok )  return(R_NilValue) ;

    //  if( nout == 0 ) return(R_NilValue);

    //  allocate nout x 3 matrix and copy from vertout to it
    SEXP    out = PROTECT( Rf_allocMatrix(REALSXP,nout,3) );
    double  *outmat = REAL(out);

    for( int i=0 ; i<nout ; i++ )
        {
        outmat[i]           = vertout[0][i] ;
        outmat[i + nout]    = vertout[1][i] ;
        outmat[i + 2*nout]  = vertout[2][i] ;
        }

    UNPROTECT(1);

    return( out );
    }
