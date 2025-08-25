
#include <R.h>
#include <Rinternals.h>

//  #include <Defn.h>               // for LDOUBLE

//  #include <float.h>
#include <stdbool.h>

//  #include "config.h"         // this file is created by configure.win or configure
//  #include "macros.h"

#include "area.h"
#include "clipping.h"

/* Required by C99, but might be slow */
#ifdef HAVE_LONG_DOUBLE
# define LDOUBLE long double
#else
# define LDOUBLE double
#endif



//  load 4 rows of quadmat
//  this is transposed from the function by the similar name in linking.c
static  bool
load_quadmat( const double offset[3], const double edge1[3], const double edge2[3],
                        double quadmat[3][4] )
    {
    for( int j=0 ; j<3 ; j++ )
        {
        quadmat[j][0]   = offset[j] - 0.5 * edge1[j] - 0.5*edge2[j] ;
        quadmat[j][1]   = offset[j] - 0.5 * edge1[j] + 0.5*edge2[j] ;
        quadmat[j][2]   = offset[j] + 0.5 * edge1[j] + 0.5*edge2[j] ;
        quadmat[j][3]   = offset[j] + 0.5 * edge1[j] - 0.5*edge2[j] ;
        }

    return( true );
    }


    
//  quadmat[][]     a quadrilateral in R^3

//  1.  clip to non-negative octant
//  2.  project onto face of octahedron (central projection), and then to x-y plane (ortho projection)
//  3.  return area of resulting polygon (if anything is left)

//  if 0 is ON the clipped quadrilateral returns NA_REAL, before any projections
//  tested situations: 0 a vertex, 0 on the edge, 0 in the interior of the quadrilateral

static  double
clip_project_measure( double quad[3][4] )
    {
    double  vertout[3][VERTSOUTMAX] ;

    double  *v[3]  = { vertout[0], vertout[1], vertout[2] } ;

    int nout;

    bool    ok;

    ok  = clipquad3D( quad, v, &nout ) ;
    if( ! ok  ||  nout == 0 )  return(0);

    //  vertout[][] now holds the clipped quad
    //  project onto face of octahedron and then x-y plane
    double  x[VERTSOUTMAX], y[VERTSOUTMAX] ;

    //  double  tol = 0;

    for( int i=0 ; i<nout ; i++ )
        {
        double  sum     = vertout[0][i] +  vertout[1][i] +  vertout[2][i] ;

        if( sum == 0 )
            //  this means that V_i is 0, so 0 is on the quadrilateral
            return(NA_REAL);

        x[i]    = vertout[0][i] / sum ;
        y[i]    = vertout[1][i] / sum ;
        }

    //  x[] and y[] define a polygon inside the triangle with these vertices:
    //      [0,0]  and  [1,0]  and  [0,1]

    return( area_polygon(x,y,nout) );
    }


//  smatgen     3xN matrix of generators, for the simplified matroid
//  sidxpair    N(N-1)/2 x 3 integer matrix of pairs of generators, 1-based
//  scenter     N(N-1)/2 x 3 matrix of pgram centers, in centered zonohedron coordinates
//  spoint      3D-point w.r.t. computing the linking number, must not be on a pgram, not verified
//                  this is in centered zonohedron coordinates
//
//  sidxpair and scenter do not have to be extended with antipodal data
//  this function does the extension automatically, whenever spoint is *not* 0 (the center of symmetry)

//  returns integer vector of length 1
//  method used: Gauss's integral definition,
//  using area of the L1 sphere in the non-negative octant only - an equilateral triangle.
//  because we only use 1 octant, symmetry cannot be used here

SEXP
linkingnumber3( SEXP smatgen, SEXP sidxpair, SEXP scenter, SEXP spoint )
    {
    const   int *dim ;

    dim = INTEGER(Rf_getAttrib(smatgen, R_DimSymbol));
    if( dim[0] != 3  ||  dim[1] < 3 )
        {
        Rprintf( "bad smatgen %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }
    int n = dim[1] ;
    const   double  *matgen = REAL(smatgen);

    int facets = (n*(n-1))/2 ;

    dim = INTEGER(Rf_getAttrib(sidxpair, R_DimSymbol));
    if( dim[0] != facets  || dim[1] != 2  )
        {
        Rprintf( "bad sidxpair %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }
    const   int *idxpair = INTEGER(sidxpair);

    dim = INTEGER(Rf_getAttrib(scenter, R_DimSymbol));
    if( dim[0] != facets  || dim[1] != 3  )
        {
        Rprintf( "bad scenter %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }
    const   double  *center = REAL(scenter);

    if( Rf_length(spoint) != 3 )    return(R_NilValue);
    const   double  *point = REAL(spoint);

    SEXP    out = PROTECT( Rf_allocVector(INTSXP,1) );
    *( INTEGER(out) ) = NA_INTEGER ;

    //  quadmat[][] holds the vertices of the pgram facet
    double  quadmat[3][4];

    double  offset[3];  // offset from point[]

    //  long double area = 0 ;   a tiny bit slower, accuracy a bit better but double is good enough
    double  area = 0 ; 
    //  Rprintf( "sizeof(area) = %d\n", sizeof(area) );
    
    for( int k=0 ; k<facets ; k++ )
        {
        offset[0]   = center[k]             - point[0] ;
        offset[1]   = center[k + facets]    - point[1];
        offset[2]   = center[k + 2*facets]  - point[2] ;

        //  i and j are 1-based
        int i = idxpair[k] ;
        int j = idxpair[k + facets] ;

        const double    *edge1  = matgen + 3*(i-1) ;    //  1-based to 0-based
        const double    *edge2  = matgen + 3*(j-1) ;    //  1-based to 0-based

        //  load the stored pgram into quadmat
        load_quadmat( offset, edge1 ,edge2, quadmat );

        double  a1 = clip_project_measure( quadmat );

        //  repeat the calculations for the antipodal pgram
        //  compute new offset, and swap edge1 and edge2
        offset[0]   = -center[k]             - point[0] ;
        offset[1]   = -center[k + facets]    - point[1] ;
        offset[2]   = -center[k + 2*facets]  - point[2] ;

        load_quadmat( offset, edge2, edge1, quadmat ) ;

        double  a2 = clip_project_measure( quadmat );

        if( a1==NA_REAL  ||  a2==NA_REAL )
            {
            Rprintf( "linkingnumber3(). INFO.  0 is ON the quadrilateral.  Returning NA.\n" );
            *( INTEGER(out) ) = NA_INTEGER ;
            UNPROTECT(1);
            return(out);
            }

        area    += a1 + a2 ;

        //  Rprintf( "k=%d   a1=%g   a2=%g\n", k, a1, a2 );
        }

    //  Rprintf( "area = %g\n", area );

    //  0.5 is the area of equilateral triangle, after projection to xy-plane
    //  change sign to be compatible with older conventions
    //  so the CIE inverted U, which is clockwise, has linknum = +1
    double  area_normalized = -area / 0.5 ;

    int     linknum =  (int) roundf( area_normalized );

    //  mathematically, area_normalized should be an integer
    //  do a tolerance check on that

    double  delta   = area_normalized - linknum;

    //  Rprintf( "linkingnumber3(). area_normalized - linknum = %e.\n", delta );

    //  this is much larger than really necessary
    //  double precision area accumulator is good enough
    double  tol = 5.e-6 ;

    if( tol < fabs(delta) )
        {
        Rprintf( "linkingnumber3(). WARN.  fabs(area_normalized - linknum(=%d)) = |%e|  >  %g (the tolerance).  Returning NA.\n",
                    linknum, delta, tol );
        linknum = NA_INTEGER ;
        }

    *( INTEGER(out) ) = linknum ;

    UNPROTECT(1);

    return(out);
    }

