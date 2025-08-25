#include <R.h>
#include <Rinternals.h>

//  #include <Defn.h>               // for LDOUBLE

//  #include <float.h>
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



//  load and unitize 4 rows of quadmat
static  bool
load_quadmat( const double offset[3], const double edge1[3], const double edge2[3],
                        double quadmat[4][3] )
    {
    int     j;

    for( j=0 ; j<3 ; j++ )
        {
        quadmat[0][j]   = offset[j] - 0.5 * edge1[j] - 0.5*edge2[j] ;
        quadmat[1][j]   = offset[j] - 0.5 * edge1[j] + 0.5*edge2[j] ;
        quadmat[2][j]   = offset[j] + 0.5 * edge1[j] + 0.5*edge2[j] ;
        quadmat[3][j]   = offset[j] + 0.5 * edge1[j] - 0.5*edge2[j] ;
        }

    //  normalize the columns of quadmat[][]
    //  project onto the unit sphere
    for( int i=0 ; i<4 ; i++ )
        {
        double  norm = 0 ;
        for( j=0 ; j<3 ; j++ )  norm += quadmat[i][j]*quadmat[i][j] ;

        if( fabs(norm) < 5.e-16 )
            return( false );

        norm    = sqrt(norm) ;

        for( j=0 ; j<3 ; j++ ) quadmat[i][j] /= norm;
        }

    return( true );
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
//  method used: Gauss's integral definition, using standard spherical area

SEXP
linkingnumber( SEXP smatgen, SEXP sidxpair, SEXP scenter, SEXP spoint )
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


    bool    symmetric = point[0]==0  &&   point[1]==0  &&  point[2]==0 ;


    //  quadmat[][] holds the vertices of the facet
    double  quadmat[4][3];

    LDOUBLE area = 0 ;

    //  Rprintf( "sizeof(area) = %d\n", sizeof(area) );

    for( int k=0 ; k<facets ; k++ )
        {
        double  offset[3];
        offset[0]   = center[k]             - point[0] ;
        offset[1]   = center[k + facets]    - point[1];
        offset[2]   = center[k + 2*facets]  - point[2] ;

        //  i and j are 1-based
        int i = idxpair[k] ;
        int j = idxpair[k + facets] ;

        const double    *edge1  = matgen + 3*(i-1) ;    //  1-based to 0-based
        const double    *edge2  = matgen + 3*(j-1) ;    //  1-based to 0-based

        if( ! load_quadmat(offset,edge1,edge2,quadmat) )
            {
            Rprintf( "linkingnumber(). The point (%g,%g,%g) (centered) is equal to a vertex of facet %d.\n",
                                point[0], point[1], point[2], k );
            Rprintf( "    The linking number is undefined; returning NA.\n" );
            UNPROTECT(1);
            return( out );
            }

        double  area_quads ;

        area_quads  = area_spherical_triangle( quadmat[1],  quadmat[3],  quadmat[0] ) +
                            area_spherical_triangle( quadmat[3],  quadmat[1],  quadmat[2] ) ;

        //  Rprintf( "k=%d  area_quad=%g\n", k, area_quads );

#if 0
        if( k == 0 )
            {
            for( j=0 ; j<3 ; j++ )
                {
                Rprintf( "%10g %10g %10g %10g\n", quadmat[0][j],  quadmat[1][j], quadmat[2][j], quadmat[3][j] ) ;
                }
            Rprintf( "area[0] = %e\n", area );
            }
#endif

        if( ! symmetric )
            {
            //  the antipodal side will NOT give the same result
            //  repeat the calculations on the antipodal side
            //  on this pass, + point instead of -
            //  edge1 and edge2 stay the same
            offset[0]   = center[k]             + point[0] ;
            offset[1]   = center[k + facets]    + point[1] ;
            offset[2]   = center[k + 2*facets]  + point[2] ;

            if( ! load_quadmat(offset,edge1,edge2,quadmat) )
                {
                Rprintf( "linkingnumber(). The point (%g,%g,%g) (centered) is equal to a vertex of pgram %d.\n",
                                    point[0], point[1], point[2], k );
                Rprintf( "    The linking number is undefined; returning NA.\n" );
                UNPROTECT(1);
                return( out );
                }

            double  area_quad = area_spherical_triangle( quadmat[1],  quadmat[3],  quadmat[0] ) +
                            area_spherical_triangle( quadmat[3],  quadmat[1],  quadmat[2] ) ;

            area_quads  += area_quad ;

            //  Rprintf( "k=%d  area_quad=%g\n", k, area_quad );
            }

        area += area_quads ;
        }

    if( symmetric )
        {
        //  the antipodal side will give the same result, so just double it
        area    *= 2.0 ;
        }

#if 0
    if( symmetric )
        {
        //  the antipodal side will give the same result, so just double it
        area *= 2.0 ;
        }
    else
        {
        //  the antipodal side will NOT give the same result
        //  repeat the calculations on the antipodal side
        for( int k=0 ; k<facets ; k++ )
            {
            double  offset[3];

            //  on this pass, + point instead of -
            offset[0]   = center[k]             + point[0] ;
            offset[1]   = center[k + facets]    + point[1] ;
            offset[2]   = center[k + 2*facets]  + point[2] ;

            //  i and j are 1-based
            int i = idxpair[k] ;
            int j = idxpair[k + facets] ;

            const double    *edge1  = matgen + 3*(i-1) ;    //  1-based to 0-based
            const double    *edge2  = matgen + 3*(j-1) ;    //  1-based to 0-based

            if( ! load_quadmat(offset,edge1,edge2,quadmat) )
                {
                Rprintf( "linkingnumber(). The point (%g,%g,%g) is equal to a vertex of facet %d.\n",
                                    point[0], point[1], point[2], k );
                UNPROTECT(1);
                return( out );
                }

            area += area_spherical_triangle( quadmat[1],  quadmat[3],  quadmat[0] ) +
                            area_spherical_triangle( quadmat[3],  quadmat[1],  quadmat[2] ) ;
            }
        }
#endif

    //  4*M_PI is the area of the sphere
    //  change sign to be compatible with older conventions
    double  area_normalized = -area / (4*M_PI) ;

    int     linknum =  (int) roundf( area_normalized );

#if 1
    double  tol = 5.e-6 ;
    if( tol < fabs(area_normalized - linknum) )
        {
        Rprintf( "linkingnumber(). WARN.  fabs(area_normalized - linknum) = |%g|  >  %g (the tolerance).  Returning NA.\n",
                    area_normalized - linknum, tol );
        linknum = NA_INTEGER ;
        }
#endif

    *( INTEGER(out) ) = linknum ;

    UNPROTECT(1);

    return(out);
    }

