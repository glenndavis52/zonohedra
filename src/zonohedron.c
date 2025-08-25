#include <stdbool.h>

#include <R.h>
#include <Rinternals.h>
#include <float.h>

#include "macros.h"



//  vec     numeric vector
//  mask    if true then skip that entry (ignore); there will only be a few of these
//  n       length of vec and mask, n>=2 not checked

bool
allequalexcept( const double vec[], const bool skip[], int n )
    {
    bool    found=false ;
    double  value=0;

    for( int k=0 ; k<n ; k++ )
        {
        if( skip[k] ) continue ;

        if( found )
            { if( vec[k] != value )   return(false); }
        else
            { found = true ; value = vec[k] ; }
        }

    return(true);
    }


//  find the largest coordinate of vec[], in absolute value
//  in case vec[] is all 0s, return a negative number
int
largestcoord( const double *vec, int n )
    {
    int     jmax = -1;
    double  absmax = 0 ;

    for( int j=0 ; j<n ; j++ )
        {
        double  abs = fabs(vec[j]);
        if( absmax < abs )
            {
            absmax  = abs ;
            jmax    = j ;
            }
        }

    return(jmax);
    }



//  shyper      a list of integer vectors, each one of them is a hyperplane of a simple matroid of rank 3.
//              so each hyperplane defines a zonogon face of a zonohedron.
//              usually there are only 2 generators spanning a parallelogram.
//
//  shypersub   an integer m-vector of 1-based indexes into shyper
//
//  scube       an n-vector giving the coordinates of the center of the 1st facet shypersub[0]
//
//  sgen        a distinguished integer that is a subset of all the hyperplanes defined by hypersub, checked
//              This generator is the edge from which the diameter is taken.
//              The diameter is from the center of this edge to the center of the antipodal edge.
//
//  sground     an integer n-vector in increasing order equal to the ground set of the matroid
//
//  snormal     an m x 3 matrix of normal vectors.  Each row is normal to the corresponding hyperplane
//              and is outward point for the facet, and serves to orient the facet
//
//  smatgen     3xn matrix of generators of the zonohedron
//
//  scrossprods n(n-1)/2 matrix as returned from allcrossproducts(), possibly with normalization
//
//  returns a list with items:
//      *)  3 x m matrix with facet centers in the columns
//      *)  3 x m matrix with vector radii in the columns
//      *)  logical m-vector with i'th value TRUE meaning that the point 0 is in the i'th facet

SEXP
beltdata( SEXP shyper, SEXP shypersub, SEXP scube, SEXP sgen, SEXP sground, SEXP snormal, SEXP smatgen, SEXP scrossprods )
    {
    int     numhypers = Rf_length(shyper);

    int     m = Rf_length(shypersub) ;
    const   int *hypersub   = INTEGER(shypersub);

    int     n = Rf_length(sground);
    const int     *ground = INTEGER(sground);

    //  make lookup table from ground points to raw index, 1..n
    int     gmax = ground[ n-1 ];
    int     *idxfromgnd = R_Calloc( gmax+1, int );

    for( int k=0 ; k<n ; k++ )
        idxfromgnd[ ground[k] ] = k+1;  // 0-based to 1-based

    //  make a copy of scube
    if( Rf_length(scube) != n )
        {
        Rprintf( "Internal Error. Rf_length(scube)=%d  !=  %d.\n", Rf_length(scube), n );
        return(R_NilValue);
        }

    //  pcube[] starts at the 0'th center, and then alternates between edge and center points
    double  *pcube = R_Calloc( n, double );
    memcpy( pcube, REAL(scube), n*sizeof(*pcube) );

    //  pcubedelta[] is the delta between consecutive center and edge points
    double  *pcubedelta = R_Calloc( n, double );

    //  hypermask[k] = true    iff the k'th generator is in the current hyperplane
    bool    *hypermask  = R_Calloc( n, bool );

    //  get distinguished generator in ground set and then convert to raw 1-based index
    int gen     = *INTEGER(sgen);
    int genidx  = idxfromgnd[ gen ];


    const   int *dim ;

    dim = INTEGER(Rf_getAttrib(snormal, R_DimSymbol));
    if( dim[0] != m ) return(R_NilValue);
    if( dim[1] != 3 ) return(R_NilValue);
    const   double  *normalmat = REAL(snormal);


    dim = INTEGER(Rf_getAttrib(smatgen, R_DimSymbol));
    if( dim[0] != 3 ) return(R_NilValue);
    if( dim[1] != n ) return(R_NilValue);
    const   double  *matgen     = REAL(smatgen);

    dim = INTEGER(Rf_getAttrib(scrossprods, R_DimSymbol));
    if( dim[0] != 3 ) return(R_NilValue);
    if( dim[1] != (n*(n-1))/2 ) return(R_NilValue);
    const   double  *crossprods = REAL(scrossprods);

    SEXP    sradmat = PROTECT( Rf_allocMatrix(REALSXP,3,m) );
    double  *radmat = REAL(sradmat);
    memset( radmat, 0, 3*m*sizeof(*radmat) ) ;

    SEXP    scentermat = PROTECT( Rf_allocMatrix(REALSXP,3,m) );
    double  *centermat = REAL(scentermat);
    memset( centermat, 0, 3*m*sizeof(*centermat) ) ;

    SEXP    sfacet0 = PROTECT( Rf_allocVector(LGLSXP,m) ) ;
    int     *facet0 = LOGICAL(sfacet0) ;
    memset( facet0, 0, m*sizeof(*facet0) ) ;


    //  compute first center = matgen %*% pcube
    double  *center  = centermat ;      //  point to first column of centermat
    for( int i=0 ; i<3 ; i++ )
        {
        //  center[i] = 0;      not necessary
        for( int k=0 ; k<n ; k++ )
            center[i] += matgen[ i + k*3 ] * pcube[k] ;
        }

    double  edge[3] = {0,0,0} ;   //  used for the center of each edge joining 2 adjacent facets

    for( int i=0 ; i<m ; i++ )
        {
        //  get the i'th hyperplane.
        int         hyperidx = hypersub[i] ;

        if( hyperidx<1  ||  numhypers<hyperidx )
            {
            Rprintf( "Internal Error. hyperidx=%d  is invalid.  numhypers=%d.\n",
                        hyperidx, numhypers );
            UNPROTECT(3);   // variables sradmat, scentermat, sfacet0                        
            return(R_NilValue);
            }

        SEXP        svec = VECTOR_ELT(shyper,hyperidx-1) ;   //  convert from 1-based to 0-based
        const int   *vec = INTEGER( svec ) ;
        int         nhyp = Rf_length( svec );   //  nhyp is the number of points in the hyperplane

        //  get the i'th normal, taken from the i'th row of normalmat[]
        double  normal[3] ;
        normal[0]   = normalmat[i] ;
        normal[1]   = normalmat[i + m];
        normal[2]   = normalmat[i + 2*m];

        //  find the largest coordinate of normal[]
        int     jmax = -1;
        double  absmax=0 ;

        for( int j=0 ; j<3 ; j++ )
            {
            double  abs = fabs(normal[j]);
            if( absmax < abs )
                {
                absmax  = abs ;
                jmax    = j ;
                }
            }

        if( absmax == 0 )
            {
            Rprintf( "Internal Error. hyperidx=%d  normal vector = 0 is invalid.\n", hyperidx );
            UNPROTECT(3);   // variables sradmat, scentermat, sfacet0            
            return(R_NilValue);
            }


        memset( hypermask, 0, n*sizeof(*hypermask) );

        //  get pointer to the i'th column of the radius matrix
        double  *radvec = radmat + i*3 ;

        bool    match = false;

        for( int k=0 ; k<nhyp ; k++ )
            {
            //  get the k'th point in the hyperplane and convert to raw 1-based index genkidx
            int     genkidx = idxfromgnd[ vec[k] ];

            hypermask[ genkidx-1 ]  = true ;    //  mark genkidx as being a generator of this hyperplane

            //  compare normal to the crossproduct
            int     diffidx = genidx - genkidx ;
            int     pairidx;

            if( 0 < diffidx )
                pairidx =   PAIRINDEX( genkidx, genidx, n ) - 1;   // 1-based to 0-based
            else if( diffidx < 0 )
                pairidx =   PAIRINDEX( genidx, genkidx, n ) - 1;   // 1-based to 0-based
            else
                {
                //  record and then ignore exact match
                match   = true;
                pcubedelta[genidx-1]    = 0 ;
                continue ;
                }

#if 1
            if( pairidx<0  ||  (n*(n-1))/2 <= pairidx )
                {
                Rprintf( "Internal Error. i=%d  k=%d  pairidx=%d.  genidx=%d  genkidx=%d.\n",
                            i, k, pairidx, genidx, genkidx );
                UNPROTECT(3);   // variables sradmat, scentermat, sfacet0                            
                return(R_NilValue);
                }
#endif

            double  signhalf;
            signhalf    = diffidx * normal[jmax] * crossprods[jmax + pairidx*3] ;

            //  the generators of this facet are LI, so signhalf cannot be 0
            signhalf    = 0<signhalf  ?  0.5  :  -0.5 ;

            //  record signhalf to be used later
            pcubedelta[genkidx-1] = signhalf ;

            //  increment radvec
            for( int j=0 ; j<3 ; j++ )
                radvec[j] += signhalf * matgen[j + (genkidx-1)*3] ;   // 1-based to 0-based
            }

        if( ! match )
            {
            Rprintf( "Internal Error. i=%d  distinguished point %d not found in hyperplane %d.\n",
                        i, gen, hyperidx );
            UNPROTECT(3);   // variables sradmat, scentermat, sfacet0                        
            return(R_NilValue);
            }

        if( 0 < i )
            {
            //  using radvec, advance from previous edge to i'th center
            center      = centermat + 3*i ;
            center[0]   = edge[0] + radvec[0];
            center[1]   = edge[1] + radvec[1];
            center[2]   = edge[2] + radvec[2];

            for( int k=0 ; k<nhyp ; k++ )
                {
                //  get the k'th point in the hyperplane and convert to raw 1-based index genkidx
                int     genkidx = idxfromgnd[ vec[k] ];

                //  advance pcube so it maps to i'th center
                pcube[genkidx-1] += pcubedelta[genkidx-1] ;
                }
            }


#if 0
        //  dump for tracing
        Rprintf( "mask:" );
        for( int k=0 ; k<n ; k++ )
            Rprintf( " %5d", hypermask[k] );
        Rprintf( "\n" );

        Rprintf( "cube:" );
        for( int k=0 ; k<n ; k++ )
            Rprintf( " %5g", pcube[k] );
        Rprintf( "\n" );
#endif

#if 0
        //  pcube now maps to the center of a facet
        //  verify that all pcube values are +-0.5 except for the generators of the facet, which are 0
        int badcount=0 ;
        for( int k=0 ; k<n ; k++ )
            {
            bool    ok = hypermask[k] ? (pcube[k]==0) : (fabs(pcube[k])==0.5) ;
            if( ! ok )  badcount++ ;
            }

        if( 0 < badcount )
            {
            Rprintf( "Internal Error. i=%d  pcube has %d bad values.\n", i, badcount );
            UNPROTECT(3);   // variables sradmat, scentermat, sfacet0            
            return(R_NilValue);
            }
#endif



        if( allequalexcept( pcube, hypermask, n ) )
            {
            facet0[ i ] = true ;    //  add current facet to facet0
            }

        //  using radvec, advance from i'th center to next edge
        edge[0]   = center[0] + radvec[0];
        edge[1]   = center[1] + radvec[1];
        edge[2]   = center[2] + radvec[2];

        for( int k=0 ; k<nhyp ; k++ )
            {
            //  get the k'th point in the hyperplane and convert to raw 1-based index genkidx
            int     genkidx = idxfromgnd[ vec[k] ];

            //  advance pcube so it maps to the next edge
            pcube[genkidx-1] += pcubedelta[genkidx-1] ;
            }
        }


    R_Free( pcube );
    R_Free( pcubedelta );
    R_Free( idxfromgnd );
    R_Free( hypermask );


    SEXP    out = PROTECT( Rf_allocVector(VECSXP,3) );
    SET_VECTOR_ELT( out, 0, sradmat );
    SET_VECTOR_ELT( out, 1, scentermat );
    SET_VECTOR_ELT( out, 2, sfacet0 );

    UNPROTECT(4);   // variables sradmat, scentermat, sfacet0, and  out

    return( out );
    }


//  shyper      a list of integer vectors, each one of them is a hyperplane of a simple matroid of rank 3.
//              so each hyperplane defines a zonogon face of a zonohedron.
//              usually there are only 2 generators spanning a parallelogram.
//              The integers are all positive and in the ground set; see the next argument sground
//
//  sground     an integer n-vector in increasing order equal to the ground set of the matroid
//
//  sradiusgen  a numeric n-vector with the radius of each generator, using raw 0-based index for the generators
//
//  returns a numeric vector with the computed radius of each hyperplane, therefore with the same length as shyper

SEXP
radiusfacet( SEXP shyper, SEXP sground, SEXP sradiusgen )
    {
    int     n = Rf_length(sground);

    if( Rf_length(sradiusgen) != n )
        {
        Rprintf( "Internal Error. Rf_length(sradiusgen)=%d  !=  %d = Rf_length(sground).\n", Rf_length(sradiusgen), n );
        return(R_NilValue);
        }


    //  make lookup table from ground points to raw index, 0..n-1
    const int     *ground = INTEGER(sground);

    int     gmax = ground[ n-1 ];
    int     *idxfromgnd = R_Calloc( gmax+1, int );

    for( int k=0 ; k<n ; k++ )
        idxfromgnd[ ground[k] ] = k;  // 0-based

    const   double  *radgen = REAL(sradiusgen);

    int     m = Rf_length(shyper);

    SEXP    out = PROTECT( Rf_allocVector(REALSXP,m) );
    double  *radfac = REAL(out);

    for( int k=0 ; k<m ; k++ )
        {
        const   SEXP    svec = VECTOR_ELT(shyper,k) ;
        const int       *hyper  = INTEGER( svec ) ;
        int             nhyper  = Rf_length( svec );

        if( nhyper == 2 )
            {
            //  this is the usual case, avoid for() loop overhead
            radfac[k]   = radgen[ idxfromgnd[hyper[0]] ]  +  radgen[ idxfromgnd[hyper[1]] ] ;
            }
        else
            {
            //  a non-trivial hyperplane, which is rare
            radfac[k]   = 0;
            for( int j=0 ; j<nhyper ; j++ )
                radfac[k]   += radgen[ idxfromgnd[hyper[j]] ] ;
            }
        }

    R_Free( idxfromgnd );

    UNPROTECT(1);   //  variable out

    return( out );
    }


//  shyper      a list of integer vectors, each one of them is a hyperplane of a simple matroid of rank 3.
//              so each hyperplane defines a zonogon face of a zonohedron.
//              Only the hyperplanes that intersect the plane <x,normal> = beta are in this list.
//              Let the length of this list be M.
//              There might be duplicates if the plane intersects both facets of an antipodal facet pair.
//              Usually there are only 2 generators spanning a parallelogram.
//              The integers are all positive and in the ground set; see the next argument sground
//
//  sfacetcenter    Mx3 matrix with the centers of the facets
//
//  sfacetnormal    Mx3 matrix with the outward-pointing normal of the facets
//                  This serves to orient the facet.
//
//  scenternormal   M-vector with the precomputed dot product of the facet centers and normal[]
//
//  sbeta           the scalar beta
//
//  sground         an integer N-vector in increasing order equal to the ground set of the matroid
//                  N is the number of generators.
//
//  sgennormal      N-vector with the precomputed dot product of the generators and normal[]
//
//  smatgen         3xN matrix with generators in the columns
//
//  scrossprods     3 x N(N-1)/2 matrix with pairwise crossproducts of the generators, possibly with normalization
//
//  return value:   Mx3 matrix where the i'th row is the intersection of the plane and the boundary polygon of the facet.
//                  There are 2 such points, and we choose the one where the transition is from neg to pos side,
//                  when the polygon is traversed counter-clockwise.

SEXP
sectionzonohedron( SEXP shyper, SEXP sfacetcenter, SEXP sfacetnormal, SEXP scenternormal, SEXP sbeta,
                        SEXP sground, SEXP sgennormal, SEXP smatgen, SEXP scrossprods )
    {
    //  verify dimensions related to m = number of facets
    int m = Rf_length(shyper) ;

    const   int *dim ;
    bool    ok;

    dim = INTEGER(Rf_getAttrib(sfacetcenter, R_DimSymbol));
    ok  = dim[0]==m  &&  dim[1]==3 ;
    if( ! ok )  return(R_NilValue) ;


    dim = INTEGER(Rf_getAttrib(sfacetnormal, R_DimSymbol));
    ok  = dim[0]==m  &&  dim[1]==3 ;
    if( ! ok )  return(R_NilValue) ;


    if( Rf_length(scenternormal) != m ) return(R_NilValue) ;


    //  verify dimensions related to n = number of generators
    int n = Rf_length(sground) ;

    if( Rf_length(sgennormal) != n )    return(R_NilValue) ;

    dim = INTEGER(Rf_getAttrib(smatgen, R_DimSymbol));
    ok  = dim[0]==3  &&  dim[1]==n ;
    if( ! ok )  return(R_NilValue) ;

    dim = INTEGER(Rf_getAttrib(scrossprods, R_DimSymbol));
    ok  = dim[0]==3  &&  dim[1]==(n*(n-1))/2 ;
    if( ! ok )  return(R_NilValue) ;

    if( Rf_length(sbeta) != 1 ) return(R_NilValue);


    //  get pointers to the vectors and matrices
    const   double  *facetcenter = REAL(sfacetcenter);
    const   double  *facetnormal = REAL(sfacetnormal);
    const   double  *gennormal   = REAL(sgennormal);
    const   double  *crossprods = REAL(scrossprods);
    const   double  *centernormal = REAL(scenternormal);
    const   double  *matgen = REAL(smatgen);

    double  beta    = *(REAL(sbeta)) ;

    //  make lookup table from ground points to raw index, 0..n-1
    const int     *ground = INTEGER(sground);

    int     gmax = ground[ n-1 ];
    int     *idxfromgnd = R_Calloc( gmax+1, int );

    for( int k=0 ; k<n ; k++ )
        idxfromgnd[ ground[k] ] = k;  // 0-based index into columns of smatgen


    //  load the coefficients that traverse a parallelogram first by generator 0 and then generator 1.
    //  we call this the "standard" order
    double  vertexcoeff[4][2] = { {-0.5,-0.5}, {0.5,-0.5}, {0.5,0.5}, {-0.5,0.5} };

    //  scratch vector that holds dot product of the parallelogram vertices with the normal
    double  value[4] ;

    const   int     next4[4] = { 1, 2, 3, 0 } ;

    SEXP    out = PROTECT( Rf_allocMatrix(REALSXP,m,3) );
    double  *outmat = REAL(out);

    for( int i=0 ; i<m*3 ; i++ )
        outmat[i]   = NA_REAL ;

    for( int k=0 ; k<m ; k++ )
        {
        SEXP        svec = VECTOR_ELT(shyper,k) ;
        int         nhyp = Rf_length( svec );   //  nhyp is the number of points in the hyperplane

        if( nhyp != 2 ) continue ;              //  ignore the non-trivial hyperplanes, they will be handled in the R calling function.

        const int   *hyper = INTEGER( svec ) ;

        //  get the facet normal from the precomputed array
        double  normalfacet[3] ;
        normalfacet[0]  = facetnormal[ k ];
        normalfacet[1]  = facetnormal[ k + m ];
        normalfacet[2]  = facetnormal[ k + 2*m ];

        int jmax = largestcoord( normalfacet, 3 );
        if( jmax < 0 )
            {
            Rprintf( "Internal Error. k=%d  normal facet vector = 0 is invalid.\n", k );
            UNPROTECT(1);   // variable out         
            return(R_NilValue);
            }

        //  get the 2 generators and orient them according to outward-pointing normalfacet[]
        int j0  = idxfromgnd[ hyper[0] ] ;
        int j1  = idxfromgnd[ hyper[1] ] ;

        if( j1 < j0 )
            {
            //  we always maintain increasing order in the hyperplanes,
            //  so this should never happen, but handle it anyway.
            //  swap j0 and j1
            int temp=j0 ; j0=j1 ; j1=temp ;
            }

        //  get the crossproduct of generators j0 and j1 from precomputed array
        int pairidx =   PAIRINDEX( j0+1, j1+1, n ) - 1;   // 1-based to 0-based
        const   double  *crossprod  = crossprods + pairidx*3 ;

        if( normalfacet[jmax] * crossprod[jmax] < 0 )
            {
            //  swap j0 and j1
            int temp=j0 ; j0=j1 ; j1=temp ;
            }

        //  we now know that traversing the parallelogram in the "standard" order is counterclockwise
        //  compute the values at the vertices in this order
        for( int i=0 ; i<4 ; i++ )
            {
            value[i] = vertexcoeff[i][0] * gennormal[j0]  +  vertexcoeff[i][1] * gennormal[j1]  +  centernormal[k] - beta ;
            }

        //  traverse value[] and check that there is exactly 1 transition from neg to non-neg
        int itrans=-1;
        int count=0 ;
        for( int i=0 ; i<4 ; i++ )
            {
            if( value[i]<0  &&  0<=value[ next4[i] ] )
                {
                itrans = i ;  count++ ;
                }
            }

        if( count != 1 )
            {
            Rprintf( "Internal Error.  k=%d.  value[] has %d transitions, but expected 1.\n", k, count );
            Rprintf( "    %g %g %g %g\n", value[0], value[1], value[2], value[3] );
            continue ;
            //return(R_NilValue);
            }

        //  go back to value[], and the actual generators, to compute the intersection on *right* side of parallelogram
        const   double  *gen0   = matgen + j0*3 ;
        const   double  *gen1   = matgen + j1*3 ;
        int     i0 = itrans ;
        int     i1 =  next4[i0] ;

        double  lambda = value[i1] / (value[i1] - value[i0]);

        for( int j=0 ; j<3 ; j++ )
            {
            double  x0 = vertexcoeff[i0][0] * gen0[j]  +  vertexcoeff[i0][1] * gen1[j] ;
            double  x1 = vertexcoeff[i1][0] * gen0[j]  +  vertexcoeff[i1][1] * gen1[j] ;

            outmat[k + j*m] = lambda*x0  +  (1-lambda)*x1  +  facetcenter[k + j*m] ;
            }
        }

    R_Free( idxfromgnd );

    UNPROTECT(1);   //  variable out

    return( out );
    }



//  shyper      a list of integer vectors, each one of them is a hyperplane of a simple matroid of rank 3.
//              so each hyperplane defines a zonogon face of a zonohedron.
//              usually there are only 2 generators spanning a parallelogram.
//
//  shypersub   an integer m-vector of 1-based indexes into shyper
//
//  sgen        a distinguished integer that is an element of all the hyperplanes defined by hypersub, checked
//              This generator is the edge from which the diameter is taken.
//              The diameter is from the center of this edge to the center of the antipodal edge.
//
//  scenter     numhypers x 3 matrix with the centers of the m facets in the rows
//
//  snormal     numhypers x 3 matrix of normal vectors.  Each row is normal to the corresponding hyperplane
//              and is outward point for the facet, and serves to orient the facet
//  sground     an integer n-vector in increasing order equal to the ground set of the matroid
//
//  smatgen     3xn matrix of generators of the zonohedron
//
//  scrossprods n(n-1)/2 matrix as returned from allcrossproducts(), possibly with normalization
//
//  returns:
//      *)  m x 3 matrix with midpoints of the sgen edges in the columns


SEXP
beltmidpoints( SEXP shyper, SEXP shypersub, SEXP sgen, SEXP scenter, SEXP snormal, SEXP sground, SEXP smatgen, SEXP scrossprods )
    {
#if 0
    Rprintf( "TYPEOF(shyper)=%d.  TYPEOF(0)=%d\n", TYPEOF(shyper), TYPEOF(VECTOR_ELT(shyper,0)) );
    Rprintf( "TYPEOF(shypersub) = %d.\n", TYPEOF(shypersub) );
    Rprintf( "TYPEOF(sgen) = %d.\n", TYPEOF(sgen) );
    Rprintf( "TYPEOF(sground) = %d.\n", TYPEOF(sground) );

    Rprintf( "TYPEOF(scenter)=%d    TYPEOF(dim)=%d.\n", TYPEOF(scenter), TYPEOF(Rf_getAttrib(scenter,R_DimSymbol)) );
    Rprintf( "TYPEOF(snormal)=%d    TYPEOF(dim)=%d.\n", TYPEOF(snormal), TYPEOF(Rf_getAttrib(snormal,R_DimSymbol)) );
    Rprintf( "TYPEOF(smatgen)=%d    TYPEOF(dim)=%d.\n", TYPEOF(smatgen), TYPEOF(Rf_getAttrib(smatgen,R_DimSymbol)) );
    Rprintf( "TYPEOF(scrossprods)=%d    TYPEOF(dim)=%d.\n", TYPEOF(scrossprods), TYPEOF(Rf_getAttrib(scrossprods,R_DimSymbol)) );
#endif

    int     numhypers = Rf_length(shyper);

    int     m = Rf_length(shypersub) ;
    const   int *hypersub   = INTEGER(shypersub);

    int     n = Rf_length(sground);
    const int     *ground = INTEGER(sground);


    //  check dimensions
    const   int *dim ;

    dim = INTEGER(Rf_getAttrib(scenter, R_DimSymbol));
    if( dim[0] != numhypers ||  dim[1] != 3 )
        {
        Rprintf( "Internal Error. center bad dimensions %dx%d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }

    dim = INTEGER(Rf_getAttrib(snormal, R_DimSymbol));
    if( dim[0] != numhypers  ||  dim[1] != 3 )
        {
        Rprintf( "Internal Error. normal bad dimensions %dx%d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }

    dim = INTEGER(Rf_getAttrib(smatgen, R_DimSymbol));
    if( dim[0] != 3  ||  dim[1] != n )
        {
        Rprintf( "Internal Error. matgen bad dimensions %dx%d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }

    dim = INTEGER(Rf_getAttrib(scrossprods, R_DimSymbol));
    if( dim[0] != 3  ||  dim[1] != (n*(n-1))/2 )
        {
        Rprintf( "Internal Error. crossprods bad dimensions %dx%d.\n", dim[0], dim[1] );
        return(R_NilValue);
        }

    //  make lookup table from ground points to raw index, 1..n
    int     gmax = ground[ n-1 ];
    int     *idxfromgnd = R_Calloc( gmax+1, int );

    for( int k=0 ; k<n ; k++ )
        idxfromgnd[ ground[k] ] = k+1;  // 0-based to 1-based


    //  get distinguished generator in ground set and then convert to raw 1-based index
    int gen     = *(INTEGER(sgen));
    int genidx  = idxfromgnd[ gen ];


    //  get pointers
    const   double  *centermat  = REAL(scenter);
    const   double  *normalmat  = REAL(snormal);
    const   double  *matgen     = REAL(smatgen);
    const   double  *crossprods = REAL(scrossprods);

    //  allocate output
    SEXP    out = PROTECT( Rf_allocMatrix(REALSXP,m,3) );
    double  *midpointmat = REAL(out);
    //  memset( midpointmat, 0, 3*m*sizeof(*midpointmat) ) ;

    for( int i=0 ; i<m ; i++ )
        {
        //  get the i'th hyperplane.
        int         hyperidx = hypersub[i] ;    //  hyperidx is 1-based

        if( hyperidx<1  ||  numhypers<hyperidx )
            {
            Rprintf( "Internal Error. hyperidx=%d  is invalid.  numhypers=%d.\n",
                        hyperidx, numhypers );
            UNPROTECT(1);   // variable out                        
            return(R_NilValue);
            }
            
        SEXP        svec = VECTOR_ELT(shyper,hyperidx-1) ;   //  convert from 1-based to 0-based
        
        //Rprintf( "i=%d   TYPEOF(shyper,%d)=%d   nhyp=%d.\n", i, hyperidx-1, TYPEOF(svec), Rf_length(svec) );          
        //continue ;
        
        const int   *vec = INTEGER( svec ) ;
        int         nhyp = Rf_length( svec );   //  nhyp is the number of points in the hyperplane
        
        //  get the i'th normal, taken from the hyperidx'th row of normalmat[]
        double  normal[3] ;
        normal[0]   = normalmat[hyperidx-1] ;
        normal[1]   = normalmat[hyperidx-1  +  numhypers];
        normal[2]   = normalmat[hyperidx-1  +  numhypers*2];

        //  find the largest coordinate of normal[]
        int jmax = largestcoord( normal, 3 );
        if( jmax < 0 )
            {
            Rprintf( "Internal Error. i=%d  normal facet vector = 0 is invalid.\n", i );
            UNPROTECT(1);   // variable out            
            return(R_NilValue);
            }

        double  radius[3] = { 0, 0, 0 };

        bool    match = false;

        for( int k=0 ; k<nhyp ; k++ )
            {
            //  get the k'th point in the hyperplane and convert to raw 1-based index genkidx
            int     genkidx = idxfromgnd[ vec[k] ];

            //  compare normal to the crossproduct
            int     diffidx = genidx - genkidx ;
            int     pairidx;

            if( 0 < diffidx )
                pairidx =   PAIRINDEX( genkidx, genidx, n ) ;   //  all indexes 1-based
            else if( diffidx < 0 )
                pairidx =   PAIRINDEX( genidx, genkidx, n ) ;   //  all indexes 1-based
            else
                {
                //  record and then ignore exact match
                match   = true;
                continue ;
                }

#if 1
            if( pairidx<1  ||  (n*(n-1))/2 < pairidx )
                {
                Rprintf( "Internal Error. i=%d  k=%d  pairidx=%d.  genidx=%d  genkidx=%d.\n",
                            i, k, pairidx, genidx, genkidx );
                UNPROTECT(1);   // variable out                            
                return(R_NilValue);
                }
#endif

            double  signhalf;
            signhalf    = diffidx * normal[jmax] * crossprods[jmax + (pairidx-1)*3] ; // 1-based to 0-based

            //  the generators of this facet are LI, so signhalf cannot be 0
            signhalf    = 0<signhalf  ?  0.5  :  -0.5 ;

            //  increment radius
            for( int j=0 ; j<3 ; j++ )
                radius[j] += signhalf * matgen[j + (genkidx-1)*3] ;   // 1-based to 0-based
            }

        if( ! match )
            {
            Rprintf( "Internal Error. i=%d  distinguished point %d not found in hyperplane %d.\n",
                        i, gen, hyperidx );
            UNPROTECT(1);   // variable out                        
            return(R_NilValue);
            }

        //  get the i'th center, taken from the hyperidx'th row of centermat[]
        double  center[3] ;
        center[0]   = centermat[hyperidx-1] ;
        center[1]   = centermat[hyperidx-1  +  numhypers];
        center[2]   = centermat[hyperidx-1  +  numhypers*2];
        
        //  midpoint is center + radius        
        midpointmat[i]         = center[0] + radius[0] ;
        midpointmat[i + m]     = center[1] + radius[1] ;
        midpointmat[i + m*2]   = center[2] + radius[2] ;
        }

    R_Free( idxfromgnd );

    UNPROTECT(1);   //  variable out

    return( out );
    }



//  this function is designed to be called many times,
//  with sdestmat and sdiff the same for each call, and modified in-place.
//  ssrcmat and sdestidx are different for each call, and read-only.
//  Each row of sdesmat is only written once,
//  and sdiff records the differences in case of *attempted* overwrites. 

//  sdestmat    N x D matrix of doubles, which is thought of as N vectors of dimension D.
//              This is the primary destination, and for the 1st call it must be initialized to NA_real_
//              This variable is both read and write;  it is modified in-place.
//
//  sdiff       real N-vector for recording the discrepancy in each row of sdestmat,
//              and for the 1st call it must be initialized to all 0s.
//              This variable is both read and write;  it is modified in-place.
//
//  ssrcmat     M x D matrix of doubles, which is viewed as M vectors of dimension D.
//              These vectors are copied to sdestmat.
//
//  sdestidx    integer M-vector with 1-based indexes into sdestmat.
//              All values must be between 1 and N.
//              Each index is to a row of destmat, and is the destination row of the corresponding
//              source row of srcmat.
//
//  returns the number of rows actually copied, or NULL in case of error
//
//  The first time this is called, it will return M meaning all rows copied.

SEXP
multicopy( SEXP sdestmat, SEXP sdiff, SEXP ssrcmat, SEXP sdestidx )
    {
    const   int *dim ;

    dim = INTEGER(Rf_getAttrib(sdestmat, R_DimSymbol));
    int n = dim[0];
    int d = dim[1];

    if( Rf_length(sdiff) != n ) return(R_NilValue);


    dim = INTEGER(Rf_getAttrib(ssrcmat, R_DimSymbol));
    int m = dim[0];
    if( dim[1] != d ) return(R_NilValue);


    if( Rf_length(sdestidx) != m ) return(R_NilValue);

    double          *destmat = REAL(sdestmat);
    double          *diff = REAL(sdiff);
    const double    *srcmat = REAL(ssrcmat);
    const int       *destidx = INTEGER(sdestidx);

    int     numcopies = 0;

    for( int i=0 ; i<m ; i++ )
        {
        int idest   = destidx[i];

#if 1
        if( idest<1 || n<idest )
            {
            Rprintf( "multicopy().  destidx[%d] = %d is invalid.\n", i, destidx[i] );
            return( R_NilValue );
            }
#endif
        double          *destvec    = destmat + idest-1 ;  // 1-based to 0-based
        const   double  *srcvec     = srcmat + i ;

        if( R_IsNA(destvec[0]) )
            {
            //  first time do the copy
            for( int j=0 ; j<d ; j++ )
                destvec[j*n]    = srcvec[j*m] ;

            numcopies++ ;
            }
        else
            {
            //  row idest has already been copied for the first time
            //  compare src and dest and record the difference
            for( int j=0 ; j<d ; j++ )
                {
                double  delta = fabs(destvec[j*n] - srcvec[j*m]) ;

                diff[idest-1]   = MAX2( delta, diff[idest-1] ) ;
                }
            }
        }

    //  return numcopies

    SEXP    out = PROTECT( Rf_allocVector(INTSXP,1) );
    *(INTEGER(out)) = numcopies;

    UNPROTECT(1);   //  variable out

    return(out);
    }