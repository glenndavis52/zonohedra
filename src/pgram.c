
#include <R.h>
#include <Rinternals.h>

#include <stdbool.h>
#include <float.h>

#include "macros.h"

//  now some more macros

#define SIGNOF(x)   ( 0<(x)  ?  1  :  ( (x)<0 ? -1 : 0 ) )

//  in this one, b must be positive
//  the interval [-b,b] maps to 0
#define SQUEEZE(x,b)    ( (x)<(-(b))  ?  (x)+(b)  :  ( (b)<(x) ? (x)-(b) : 0 ) )

#define DOT2(v,w)   ( v[0]*w[0]  +  v[1]*w[1] )




//  p1, p2  points in the plane

//  returns distance between the points squared
static  inline  double
distbetween2( const double p1[2], const double p2[2] )
    {
    double  diff[2] = { p1[0] - p2[0], p1[1] - p2[1] } ;

    return( DOT2(diff,diff) );
    }




//  p   the given point in the plane
//  g1  1st generator of pgram
//  g2  2nd generator of pgram, these must be linearly independent

//  returns distance^2 from p to the pgram

double
dist2pgramSQ_2D( const double p[2], const double g1[2], const double g2[2] )
    {
    //  copy generators to local variables
    double  v1[2] = {g1[0],g1[1]} ;
    double  v2[2] = {g2[0],g2[1]} ;

    double  v1v2    = DOT2(v1,v2) ;

    if( v1v2 < 0 )
        {
        //  force acute angle between v1 and v2, by changing sign of v2
        v2[0]   *= -1 ;
        v2[1]   *= -1 ;
        v1v2    *= -1 ;
        }

    //double  v2v1 = v1v2 ;
    //double  v1v1 = DOT2(v1,v1) ;
    //double  v2v2 = DOT2(v2,v2) ;


    double  len ;

    //  unit normals inward pointing

    double  u1[2] = { v1[1], -v1[0] } ;
    len     = sqrt( DOT2(u1,u1) ) ;
    u1[0]   /= len ;
    u1[1]   /= len ;
    double  u1v2 = DOT2(u1,v2) ;
    if( u1v2 < 0 )
        {
        u1[0]   *= -1 ;
        u1[1]   *= -1 ;
        u1v2    *= -1 ;
        }


    double  u2[2] = { v2[1], -v2[0] } ;
    len     = sqrt( DOT2(u2,u2) ) ;
    u2[0]   /= len ;
    u2[1]   /= len ;
    double  u2v1 = DOT2(u2,v1) ;
    if( u2v1 < 0 )
        {
        u2[0]   *= -1 ;
        u2[1]   *= -1 ;
        u2v1    *= -1 ;
        }

    //  compute limits of the slabs
    double  b1 = DOT2(u1,v2) / 2 ;
    double  b2 = DOT2(u2,v1) / 2 ;

    //  compute parameters of the input point p[]
    double  pu1 = DOT2(p,u1) ;
    double  pu2 = DOT2(p,u2) ;

    //  squeeze both
    double  t1  = SQUEEZE(pu1,b1) ;
    double  t2  = SQUEEZE(pu2,b2) ;

    //  compute the signs, +1 or -1 or 0
    int s1 = SIGNOF(t1) ;
    int s2 = SIGNOF(t2) ;

    //  from the 2 signs compute the region, and index from -4 to 4
    int region = s1 + 3*s2 ;

    if( region == 0 )
        //  both signs are 0
        //  p is inside the pgram, the trivial case
        return( 0 );

    //  precompute some dot products, that are potentially needed
    double  pv1 = DOT2(p,v1);
    double  pv2 = DOT2(p,v2);

    //  compute vertex
    //  this might be used because
    //      *) it is the closest point to p
    //      *) it is used to compute the boundary of a "fan" for a decision threshold

    double  vertex[2] ;

    //  compute special signs for the vertex
    //  NB:  these are swapped from s1 and s2
    int sv1 = s2 ;
    int sv2 = s1 ;

    //  since region!=0, both of these cannot be 0
    if( sv1 == 0 )  sv1 = sv2 ;

    if( sv2 == 0 )  sv2 = sv1 ;

    vertex[0]   = 0.5*( sv1*v1[0] + sv2*v2[0] ) ;
    vertex[1]   = 0.5*( sv1*v1[1] + sv2*v2[1] ) ;

    double  thresh1, thresh2 ;

    double  d2=0 ;

    switch( region )
        {
        case -4 :
        case  4 :
            //  both signs are the same, and the angle between v1 and v2 is acute
            //  so the closest point is a vertex
            d2  = distbetween2(p,vertex) ;
            break;


        case -1 :
            //  s1=-1 and s2=0, 1 test needed
            thresh1  = DOT2(vertex,v1) ;
            if( thresh1 < pv1 )
                //  distance2 to an edge
                d2 = t1*t1 ;
            else
                //  distance2 to a vertex
                d2  = distbetween2(p,vertex) ;
            break;

        case +1 :
            //  s1=+1 and s2=0, 1 test needed
            thresh1  = DOT2(vertex,v1) ;
            if( pv1 < thresh1 )
                //  distance2 to an edge
                d2 = t1*t1 ;
            else
                //  distance2 to a vertex
                d2  = distbetween2(p,vertex) ;
            break;


        case -3 :
            //  s1=0 and s2=-1, 1 test needed
            thresh2  = DOT2(vertex,v2) ;
            if( thresh2 < pv2 )
                //  distance2 to an edge
                d2 = t2*t2 ;
            else
                //  distance2 to a vertex
                d2  = distbetween2(p,vertex) ;
            break;

        case 3 :
            //  s1=0 and s2=+1, 1 test needed
            thresh2  = DOT2(vertex,v2) ;
            if( pv2 < thresh2 )
                //  distance2 to an edge
                d2 = t2*t2 ;
            else
                //  distance2 to a vertex
                d2  = distbetween2(p,vertex) ;
            break;

        case 2 :
            //  s1=-1 and s2=+1, 2 tests needed
            thresh1  = DOT2(vertex,v1) ;
            thresh2  = DOT2(vertex,v2) ;

            if( pv1 < thresh1 )
                //  to an edge
                d2 = t1*t1 ;
            else if( thresh2 < pv2 )
                //  to an edge
                d2 = t2*t2 ;
            else
                //  to the vertex
                d2  = distbetween2(p,vertex) ;
            break ;

        case -2 :
            //  s1=+1 and s2=-1, 2 tests needed
            thresh1  = DOT2(vertex,v1) ;
            thresh2  = DOT2(vertex,v2) ;

            if( pv2 < thresh2 )
                //  to an edge
                d2 = t2*t2 ;
            else if( thresh1 < pv1 )
                d2 = t1*t1 ;
            else
                //  to the vertex
                d2  = distbetween2(p,vertex) ;
            break;
        }

    return( d2 );
    }





//  u[3]    a unit vector
//
//  returns 3x3 rotation matrix in mat[][]
bool
rotation2pole( const double u[3], double mat[3][3] )
    {
    double  z   = u[2]<0 ? -1 : 1 ;

    //  uv  = u + pole
    double  uv[3] ;
    uv[0]   = u[0] ;
    uv[1]   = u[1] ;
    uv[2]   = u[2] + z ;

    double  s = 1 / (1 + u[2]*z) ;

    for( int i=0 ; i<3 ; i++ )
        {
        for( int j=0 ; j<3 ; j++ )
            mat[i][j]   = -s * uv[i] * uv[j] ;
        }

    for( int j=0 ; j<3 ; j++ )
        {
        mat[j][j] += 1 ;    //  add 1 to diagonal

        mat[2][j] += 2 * z * u[j] ;
        }

    return( true );
    }


//  this function is only used for testing
SEXP
rotation2pole_test( SEXP su )
    {
    const double    *u = REAL(su) ;

    double  mat[3][3] ;

    bool    ok = rotation2pole( u, mat ) ;
    if( ! ok )  return( R_NilValue );

    SEXP    out = PROTECT( Rf_allocMatrix(REALSXP,3,3) );
    double  *matout = REAL(out);

    //  now transpose
    for( int i=0 ; i<3 ; i++ )
        for( int j=0 ; j<3 ; j++ )
            matout[ i + 3*j ]   = mat[i][j] ;

    UNPROTECT(1) ;      // out

    return( out );
    }

//  point   the given point
//  v1      1st generator of pgram
//  v2      2nd generator of pgram, these must be linearly independent
//  center  center of the pgram
//  normal  unit normal of the pgram

//  returns distance from p to the pgram

double
dist2pgram( const double point[3], const double v1[3], const double v2[3], const double center[3], const double normal[3] )
    {
    double  mat[3][3] ;

    bool    ok = rotation2pole( normal, mat ) ;
    if( ! ok )  return( R_NaReal );

    //  compute point relative to the center
    double  p_centered[3] ;
    for( int k=0 ; k<3 ; k++ )  p_centered[k] = point[k] - center[k] ;

    //  rotate p_centered
    double  p_rot[3] ;
    for( int i=0 ; i<3 ; i++ )
        {
        p_rot[i] = 0 ;
        for( int j=0 ; j<3 ; j++ )
            p_rot[i] += mat[i][j] * p_centered[j] ;
        }

    //  for v1 and v2, we only need the first 2 coords, the 3rd is 0 or very close to it
    double v1_rot[2], v2_rot[2] ;
    for( int i=0 ; i<2 ; i++ )
        {
        v1_rot[i] = 0 ;
        v2_rot[i] = 0 ;
        for( int j=0 ; j<3 ; j++ )
            {
            v1_rot[i] += mat[i][j] * v1[j] ;
            v2_rot[i] += mat[i][j] * v2[j] ;
            }
        }

    double  dist_to_planeSQ  = p_rot[2]*p_rot[2] ;

    double  dist_within_planeSQ  = dist2pgramSQ_2D( p_rot, v1_rot, v2_rot ) ;

    double  dist    = sqrt( dist_to_planeSQ  +  dist_within_planeSQ );
    
    //Rprintf( "p=(%g,%g,%g)  center=(%g,%g,%g)  normal=(%g,%g,%g)  d2p=%g  dwithp=%g\n",
    //            point[0], point[1], point[2],  center[0], center[1], center[2],  normal[0], normal[1], normal[2],
    //            sqrt(dist_to_planeSQ), sqrt(dist_within_planeSQ) ) ;

    return( dist );
    }




//  smatgen     3xN matrix of generators, for the simplified matroid
//  sidxpair    N(N-1)/2 x 3 integer matrix of pairs of generators, 1-based
//  scenter     N(N-1)/2 x 3 matrix of pgram centers, in centered zonohedron coordinates
//  snormal     N(N-1)/2 x 3 matrix of unit normals to the pgrams
//  spoint      single 3D-point w.r.t. computing the linking number, must not be on a pgram, not verified
//                  this is in centered zonohedron coordinates
//
//  sidxpair and scenter do not have to be extended with antipodal data
//  this function does the extension automatically, whenever spoint is *not* 0 (the center of symmetry)

//  returns real vector of length 1

SEXP
dist2surface( SEXP smatgen, SEXP sidxpair, SEXP scenter, SEXP snormal, SEXP spoint )
    {
    const   int *dim ;

    dim = INTEGER(Rf_getAttrib(smatgen, R_DimSymbol));
    if( dim[0] != 3  ||  dim[1] < 3 )
        {
        Rprintf( "dist2surface().  bad smatgen %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }
    int n = dim[1] ;
    const   double  *matgen = REAL(smatgen);

    int facets = (n*(n-1))/2 ;

    dim = INTEGER(Rf_getAttrib(sidxpair, R_DimSymbol));
    if( dim[0] != facets  || dim[1] != 2  )
        {
        Rprintf( "dist2surface().  bad sidxpair %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }
    const   int *idxpair = INTEGER(sidxpair);

    dim = INTEGER(Rf_getAttrib(scenter, R_DimSymbol));
    if( dim[0] != facets  || dim[1] != 3  )
        {
        Rprintf( "dist2surface().  bad scenter %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }
    const   double  *center = REAL(scenter);

    dim = INTEGER(Rf_getAttrib(snormal, R_DimSymbol));
    if( dim[0] != facets  || dim[1] != 3  )
        {
        Rprintf( "dist2surface().  bad snormal %d x %d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }
    const   double  *normal = REAL(snormal);


    if( Rf_length(spoint) != 3 )
        {
        Rprintf( "dist2surface().  bad spoint length=%d.\n", Rf_length(spoint)  );
        return(R_NilValue);
        }
    const   double  *point = REAL(spoint);


    bool    symmetric = point[0]==0  &&   point[1]==0  &&  point[2]==0 ;

    double  dist_min = FLT_MAX ;

    for( int k=0 ; k<facets ; k++ )
        {
        //  i and j are 1-based
        int i = idxpair[k] ;
        int j = idxpair[k + facets] ;

        const double    *v1  = matgen + 3*(i-1) ;    //  1-based to 0-based
        const double    *v2  = matgen + 3*(j-1) ;    //  1-based to 0-based

        double  center_facet[3] ;
        center_facet[0]   = center[k] ;
        center_facet[1]   = center[k + facets] ;
        center_facet[2]   = center[k + 2*facets] ;

        double  normal_facet[3] ;
        normal_facet[0]   = normal[k] ;
        normal_facet[1]   = normal[k + facets] ;
        normal_facet[2]   = normal[k + 2*facets] ;

        double  dist    = dist2pgram( point, v1, v2, center_facet, normal_facet ) ;

        dist_min    = MIN2( dist, dist_min );

        if( dist_min == 0 )
            //  can't get any smaller than 0
            break ;

        if( ! symmetric )
            {
            //  the antipodal side will NOT give the same result
            //  repeat the calculations on the antipodal side
            //  reverse the center
            //  But normal can stay the same, and rotation2pole() will take care of the reversal !
            center_facet[0] *= -1 ;
            center_facet[1] *= -1 ;
            center_facet[2] *= -1 ;

            dist    = dist2pgram( point, v1, v2, center_facet, normal_facet ) ;

            dist_min    = MIN2( dist, dist_min );

            if( dist_min == 0 )
                //  can't get any smaller than 0
                break ;
            }

        }



    SEXP    out = PROTECT( Rf_allocVector(REALSXP,1) );
    *( REAL(out) ) = dist_min ;

    UNPROTECT(1);

    return(out);
    }

    

//  this function is only used for testing

SEXP
dist2pgram_test( SEXP spoint, SEXP sv1, SEXP sv2, SEXP scenter, SEXP snormal )
    {
    const   double  *point  = REAL(spoint);
    const   double  *v1     = REAL(sv1);
    const   double  *v2     = REAL(sv2);
    const   double  *center = REAL(scenter);
    const   double  *normal = REAL(snormal);

    double  d = dist2pgram( point, v1, v2, center, normal ) ;

    SEXP    out = PROTECT( Rf_allocVector(REALSXP,1) );

    REAL(out)[0]    = d ;

    UNPROTECT(1);

    return( out );
    }

    