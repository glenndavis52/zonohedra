#include <R.h>
#include <Rinternals.h>

//  #include <Defn.h>               // for LDOUBLE

#include <float.h>
#include <stdbool.h>

//  #include "config.h"         // this file is created by configure.win or configure

#include "macros.h"
#include "area.h"


/* Required by C99, but might be slow */
#ifdef HAVE_LONG_DOUBLE
# define LDOUBLE long double
#else
# define LDOUBLE double
#endif



//  smatcum     3 x N+1 matrix of the cumulative sum of the N generators, for the simplified matroid
//                  the first column is all 0s
//  spoint      3D-point w.r.t. computing the linking number, must not be on a pgram. not verified
//                  spoint is in *centered* zonohedron coordinates
//
//  returns integer vector of length 1
//  method used: Gauss's integral definition, using standard spherical area
//
//  compared with linkingnumber(), it has modified arguments.  It needs more work.

SEXP
linkingnumber2( SEXP smatcum, SEXP spoint )
    {
    const   int *dim ;

    dim = INTEGER(Rf_getAttrib(smatcum, R_DimSymbol));
    if( dim[0] != 3  ||  dim[1] <= 3 )
        {
        Rprintf( "bad smatcum %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }
    int n = dim[1] - 1;
    const   double  *matcum = REAL(smatcum);

    //  verify that the first column is all 0s
    bool    ok = matcum[0]==0  &&  matcum[1]==0  &&  matcum[2]==0 ;
    if( ! ok )
        {
        Rprintf( "matcum is invalid; 1st column must be 0.\n" );
        return(R_NilValue);
        }

    if( Rf_length(spoint) != 3 )    return(R_NilValue);
    const   double  *point = REAL(spoint);

    const double    *white = matcum + 3*n ;

    double  center[3] ;
    for( int k=0 ; k<3 ; k++ )  center[k] = white[k]/2 ;

    //  allocate (N+1) x 2N x 3     3D array
    double  *vertex = R_Calloc( (n+1) * (2*n) * 3, double );

    int     rowstep = 3 ;
    int     colstep = 3*(n+1) ;

    //  fill vertex with the *centered* vertices.
    //  There are N(N-1) + 2 vertices in the surface, and N(N-1) + 2N cells in the array.

    //  start with "black" and "white", which are duplicated to many cells in the array
    //  this section fills 2N cells
    for( int jj=0 ; jj<n ; jj++ )
        {
        int     jb = 2*jj ;
        int     jw = jb + (n % 2) ;     // when n is odd, there is a shift

        double  *v0 = vertex  +  0*rowstep  +  jb*colstep ;
        double  *v1 = vertex  +  n*rowstep  +  jw*colstep ;

        for( int k=0 ; k<3 ; k++ )
            {
            v0[k]   = -center[k];            //  *black* is -center
            v1[k]   =  center[k];            //  *white* is +center
            }
        }

    int row, col;

    //  this section fills N(N-1) cells with *centered* vertices
    for( int j1=0 ; j1<n-1 ; j1++ )
        {
        const double    *s1 = matcum + j1*3;

        for( int j2=j1+1 ; j2<n ; j2++ )
            {
            const double    *s2 = matcum + j2*3;

            row = j2 - j1 ;
            col = j1 + j2 ;
            double  *v = vertex  +  row*rowstep  +  col*colstep ;

            //  now the antipodal cell
            row = n - row ;
            col = (col + n) % (2*n) ;
            double  *va = vertex  +  row*rowstep  +  col*colstep ;

            for( int k=0 ; k<3 ; k++ )
                {
                v[k]    = s2[k] - s1[k] - center[k] ;
                va[k]   = -v[k] ;
                }
            }
        }


    SEXP    out = PROTECT( Rf_allocVector(INTSXP,1) );
    *( INTEGER(out) ) = NA_INTEGER ;

    bool    symmetric = point[0]==0  &&  point[1]==0  &&  point[2]==0 ;

    //  in the next pass, subtract off point[] and unitize
    //  there are N(N-1) + 2N  cells

    int     rowmax = symmetric  ?  n/2+1 : n ;      // this opimization does not decrease time much - Jul 26 2022

    for( row=0 ; row<=rowmax ; row++ )
        {
        for( col=(row%2) ; col<2*n ; col+=2 )
            {
            double  *v = vertex  +  row*rowstep  +  col*colstep ;

            double  norm=0 ;
            for( int k=0 ; k<3 ; k++ )
                {
                v[k]    -= point[k] ;
                norm    += v[k]*v[k] ;
                }

            if( fabs(norm) < 5.e-16 )
                {
                Rprintf( "linkingnumber2(). The point (%g,%g,%g) is equal to a vertex of the surface.\n",
                                    point[0], point[1], point[2] );
                R_Free( vertex );
                UNPROTECT(1);
                return( out );
                }

            norm    = sqrt(norm);

            for( int k=0 ; k<3 ; k++ )  v[k] /= norm ;
            }
        }


    //  iterate over the pgrams and compute area
    //  there are N(N-1) pgrams, unless symmetric when there are N(N-1)/2

    rowmax = symmetric  ?  n/2 : n-1 ;  // this opimization is effective - Jul 25 2022

    int     pgrams = 0 ;
    double  area = 0 ;
    for( row=1 ; row<=rowmax ; row++ )
        {
        int collim ;

        if( row<rowmax  ||  ! symmetric  ||  n%2 )
            collim  = 2*n ;
        else
            collim  = n ;

        for( col=((row+1)%2) ; col<collim ; col+=2 )
            {
            int colneg  = (col-1+2*n) % (2*n);
            int colpos  = (col+1) % (2*n);

            const double  *q0 = vertex  +  (row-1)*rowstep    +  col * colstep ;
            const double  *q1 = vertex  +  row*rowstep        +  colpos * colstep ;
            const double  *q2 = vertex  +  (row+1)*rowstep    +  col * colstep ;
            const double  *q3 = vertex  +  row*rowstep        +  colneg * colstep ;

            area += area_spherical_triangle( q1, q3, q0 )  +  area_spherical_triangle( q3, q1, q2 ) ;

            pgrams++ ;
            }
        }

    R_Free( vertex );

    if( symmetric )
        //  the antipodal side will give the same area sum, so just double it
        area    *= 2.0 ;


    //  4*M_PI is the area of the sphere
    //  change sign to be compatible with older conventions
    double  area_normalized = -area / (4*M_PI) ;

    int     linknum =  (int) roundf( area_normalized );

    *( INTEGER(out) ) = linknum ;


    //  Rprintf( "area_normalized=%g  area_normalized-linknum = %g\n", area_normalized, area_normalized-linknum );

#if 1
    int pgrams_corr = symmetric ? (n*(n-1))/2 : n*(n-1) ;

    if( pgrams != pgrams_corr )
        Rprintf( "ERROR. pgrams = %d  !=  %d (the correct value).\n", pgrams, pgrams_corr );

    double  tol = 5.e-7 ;
    if( tol < fabs(area_normalized - linknum) )
        Rprintf( "WARN. area_normalized - linknum = %g  >  %g\n", area_normalized - linknum, tol );
#endif

    UNPROTECT(1);       //  out

    return( out );
    }
