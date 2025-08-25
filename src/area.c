
#include <math.h>
//  #include <float.h>
#include <stdbool.h>

//  #include "config.h"         // this file is created by configure.win or configure

#include "macros.h"


//  compute 3x3 determinant
double
det3x3( const double a[3], const double b[3], const double c[3] )
    {
    double  out = 0 ;

    out +=  a[0] * (b[1]*c[2] - b[2]*c[1]) ;
    out -=  b[0] * (a[1]*c[2] - a[2]*c[1]) ;
    out +=  c[0] * (a[1]*b[2] - a[2]*b[1]) ;

    return( out );
    }



//  a, b, c      points on the unit 2-sphere
//  for best numerics, the shortest side should be opposite c
//
//  returns *signed* area of the triangle.
//  positive area is counterclockwise, when viewed from the outside of the sphere

double
area_spherical_triangle( const double a[3], const double b[3], const double c[3] )
    {
    double  det = det3x3( a, b, c );

    //  Rprintf( "det = %g\n", det );

    if( fabs(det) < 5.e-16 )    return( 0.0 );  //  det too small

    double  sinhalfa    = 0;
    double  sinhalfb    = 0;
    double  sinhalfc    = 0;

    for( int k=0 ; k<3 ; k++ )
        {
        double  diffa, diffb, diffc;

        diffa   = b[k] - c[k] ;
        diffb   = a[k] - c[k] ;
        diffc   = a[k] - b[k] ;

        sinhalfa    += diffa*diffa ;
        sinhalfb    += diffb*diffb ;
        sinhalfc    += diffc*diffc ;
        }

    //  the next test is not needed; it is taken care of by the determinant test above
    //  if( sinhalfa==0 || sinhalfb==0 || sinhalfc==0 ) return(0.0);

    sinhalfa    = 0.5 * sqrt( sinhalfa );
    sinhalfb    = 0.5 * sqrt( sinhalfb );
    sinhalfc    = 0.5 * sqrt( sinhalfc );

    sinhalfa    = MIN2( sinhalfa, 1 );      //  necessity unknown, no known examples
    sinhalfb    = MIN2( sinhalfb, 1 );      //  necessity unknown, no known examples

    double  coshalfa    = sqrt( 1.0 - sinhalfa*sinhalfa );
    double  coshalfb    = sqrt( 1.0 - sinhalfb*sinhalfb );

    double  sinprod = sinhalfa * sinhalfb ;
    double  cosprod = coshalfa * coshalfb ;

    double  T2  = (sinhalfa*sinhalfa + sinhalfb*sinhalfb - 2*sinprod*sinprod - sinhalfc*sinhalfc) / (2*cosprod);

    double  cosC    = T2 / sinprod ;
    cosC            = MIN2( cosC, 1 );          //  necessary, testing uncovered an example
    double  sinC    = sqrt( 1.0 - cosC*cosC );

    double  E ;     //  spherical Excess

#if 0
    //  this section uses atan2()
    //  in the next line, the x part can be +, -, or 0.
    //  but the y part is always positive, and so E is between 0 and 2*M_PI, as it should be !
    E   = 2.0 * atan2( sinprod * sinC, cosprod + T2 ) ;        //  atan2f() is the same speed !
#else
    //  this section uses atan()
    //  it has only none virtue over atan2() - it's a bit faster
    double  denom   = cosprod + T2;

    double  Ep ;
    if( denom != 0 )
        {
        Ep = 2.0 * atan( (sinprod * sinC) / denom ) ;

        if( Ep < 0 ) Ep += 2*M_PI;     // get a positive value using periodicity

        //if( 1.e-6 < fabs(E-Ep) )
        //    Rprintf( "y=%g  x=%g  E=%g   Ep=%g   E-Ep=%g\n", sinprod * sinC, cosprod + T2, E, Ep, E-Ep );
        }
    else
        {
        Ep = M_PI ; //  very special case
        }

    E   = Ep ;
#endif


#if 0
    if( E <= 0 )
        {
        Rprintf( "a=%g,%g,%g   b=%g,%g,%g   c=%g,%g,%g\n", a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2] );
        Rprintf( "sinhalfa=%g  sinhalfb=%g  sinhalfc=%g\n", sinhalfa, sinhalfb, sinhalfc ) ;
        Rprintf( "E=%g < 0. sinprod=%g  cosC=%g   sinC=%g  cosprod=%g   T2=%g   temp=%g\n",
                        E, sinprod, cosC, sinC, cosprod, T2, temp ) ;
        }
#endif

    //  we want the *signed* area
    //  E   *= (0 < det) ? 1.0 : -1.0 ;
    if( det < 0 )   E = -E ;

    //  Rprintf( "signed area = %g\n", E );

    return( E ) ;
    }


//  compute signed area of a planar polygon
//  counter-clockwise is positive

double
area_polygon( const double x[], const double y[], int n )
    {
    if( n < 3 ) return(0) ;     //  degenerate case
        
    double  area = 0 ;

    for( int i=0 ; i<n ; i++ )
        {
        int     i1 = (i+1) % n ;

        area += (x[i] - x[i1]) * (y[i1] + y[i]) ;
        }

    return( 0.5 * area );
    }


/*************      for testing only    ***********************/


#include <R.h>
#include <Rinternals.h>

//  #include <Defn.h>               // for LDOUBLE


SEXP
area_sphtri( SEXP sa, SEXP sb, SEXP sc )
    {
    SEXP    out = PROTECT( Rf_allocVector(REALSXP,1) );

    *(REAL(out)) = area_spherical_triangle( REAL(sa), REAL(sb), REAL(sc) ) ;

    UNPROTECT(1);

    return(out);
    }

