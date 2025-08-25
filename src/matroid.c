#include <stdbool.h>

//  #include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>

#include "matdat.h"
#include "macros.h"


extern  bool    collapseGroups1D( double *vec, int n, double eps ) ;



SEXP
collapseGroups1D_R( SEXP x, SEXP eps )
    {
    int     n =  Rf_length(x);
	SEXP    out = PROTECT( Rf_allocVector(LGLSXP,1) );

    *(LOGICAL(out)) = collapseGroups1D( REAL(x), n, *(REAL(eps)) ) ;

    UNPROTECT(1);

    return(out);
    }


//  sx      pointer to matrix of doubles,  prevalidated
//  seps    pointer to a small double, prevalidated
//  smargin pointer to an integer, 1 or 2, prevalidated
//
//  changes the matrix sx in-place

//  returns true or false
SEXP
conditionalAntipodal( SEXP sx, SEXP seps, SEXP smargin )
    {
    double  eps = *(REAL(seps)) ;
    //  int margin  = *(INTEGER(smargin)) ;

    //  int     eltStep ;   // step between elements of x
    //  int     vecStep ;   // step between vectors of x
    //  int     vecLen ;    // length of each vector
    //  int     nVec ;      // number of vectors

    //  Rprintf( "eps=%g  margin=%d\n", eps, margin );


	SEXP    out = PROTECT( Rf_allocVector(LGLSXP,1) );

    matdat  md = extractmatdat( sx, smargin );
    if( md.mat == NULL )
        {
        *(LOGICAL(out)) = false ;
        UNPROTECT(1);
        return(out);
        }


    for( int j=0 ; j<md.nVec ; j++ )
        {
        double  *vec = md.mat + j*md.vecStep ;

        for( int i=0 ; i<md.vecLen ; i++ )
            {
            double  v = vec[i*md.eltStep] ;

            //  Rprintf( "j=%d  i=%d  v=%g  fabs(v)=%g\n", j, i, v, fabs(v) );

            //  I do not know where fabs() is defined, but this seems to work !?
            if( fabs(v) <= eps ) continue;  //  so close to 0 that we can ignore it

            //  v is the first "significant" value in the vector
            if( v < 0 )
                {
                //  reverse sign of all of vec[]
                //  Rprintf( "reverse!\n"  );
                for( int k=0 ; k<md.vecLen ; k++ ) { vec[k*md.eltStep] = -vec[k*md.eltStep] ; }
                }
            break;
            }
        }

    *(LOGICAL(out)) = true ;

    UNPROTECT(1);

    return(out);
    }



//  sx      pointer to matrix of doubles,  prevalidated
//  smargin pointer to an integer, 1 or 2, prevalidated
//
//  changes the matrix sx in-place

//  returns true or false
SEXP
normalizeMatrix( SEXP sx, SEXP smargin )
    {
    SEXP    out = PROTECT( Rf_allocVector(LGLSXP,1) );

    matdat  md = extractmatdat( sx, smargin );
    if( md.mat == NULL )
        {
        *(LOGICAL(out)) = false ;
        UNPROTECT(1);
        return(out);
        }

    for( int j=0 ; j<md.nVec ; j++ )
        {
        double  *vec = md.mat + j*md.vecStep ;

        double  accum=0;

        for( int i=0 ; i<md.vecLen ; i++ )
            {
            double  v = vec[i*md.eltStep] ;
            accum   += v*v;
            }

        if( 0 < accum )
            {
            double  norm = sqrt(accum);

            for( int k=0 ; k<md.vecLen ; k++ ) { vec[k*md.eltStep] /= norm ; }
            }
        }

    *(LOGICAL(out)) = true ;

    UNPROTECT(1);   //  variable out

    return(out);
    }



//  sx      pointer to 3xn matrix of doubles,  prevalidated

//  returns a 3 x n(n-1)/2 matrix, with all crossproduct pairs, in order matching combinations(n,2)
//  the order of these crossproducts matches allpairs().
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
allcrossproducts( SEXP sx )
    {
    int *dim    = INTEGER(Rf_getAttrib(sx, R_DimSymbol));
    int nrow    = dim[0] ;
    int n       = dim[1] ;

    double  *mat     = REAL(sx) ;    //  pointer to the matrix of doubles

    if( mat == NULL  ||  nrow!=3 )
        {
        return(R_NilValue);
        }

    int colsout = ( n*(n-1) )/2 ;

    SEXP    out = PROTECT( Rf_allocMatrix(REALSXP,3,colsout) );
    double  *matout = REAL(out);

    int     count   = 0;
    for( int j1=0 ; j1<n-1 ; j1++ )
        {
        double  *v1 = mat + 3*j1;

        for( int j2=j1+1 ; j2<n ; j2++, count++ )
            {
            double  *v2 = mat + 3*j2;
            double  *cp = matout + 3*count ;

            cp[0] =  v1[1]*v2[2] - v1[2]*v2[1] ;
            cp[1] = -v1[0]*v2[2] + v1[2]*v2[0] ;
            cp[2] =  v1[0]*v2[1] - v1[1]*v2[0] ;
            }
        }

    UNPROTECT(1);   // variable out

    return( out );
    }




//  scrossprods     3 x N(N-1) matrix with unitized crossproduct pairs in the columns
//  shyperplane     list with strictly increasing integer vectors, defining non-trivial hyperplanes from the ground set, length M
//  scrossprodsref  3 x M matrix, with the "reference" crossproducts in the columns, one for each hyperplane
//  sground         N vector of increasing integers forming the ground set
//
//  scrossprods is modified in-place !

SEXP
snapcrossprods( SEXP scrossprods, SEXP shyperplane, SEXP scrossprodsref, SEXP sground )
    {
    int n   = Rf_length(sground) ;
    int m   = Rf_length(shyperplane);

    const   int *dim    = INTEGER(Rf_getAttrib(scrossprods,R_DimSymbol));

    if( dim[0] != 3  ||  dim[1] != n*(n-1)/2 )   return(R_NilValue);

    dim    = INTEGER(Rf_getAttrib(scrossprodsref,R_DimSymbol));

    if( dim[0] != 3  ||  dim[1] != m ) return(R_NilValue);


    //  make lookup table from ground points to raw index, 1..n
    const int     *ground = INTEGER(sground);

    int     gmax = ground[ n-1 ];
    int     *idxfromgnd = R_Calloc( gmax+1, int );

    for( int k=0 ; k<n ; k++ )
        idxfromgnd[ ground[k] ] = k+1 ;     //  1-based output

    double  *crossprods = REAL(scrossprods) ;

    const   double *crossprodsref = REAL(scrossprodsref) ;

    for( int k=0 ; k<m ; k++ )
        {
        const   double *cpref = crossprodsref + 3*k ;

        //  find the largest coordinate of cpref[]
        int     imax = -1;
        double  absmax=0 ;

        for( int i=0 ; i<3 ; i++ )
            {
            double  abs = fabs(cpref[i]);
            if( absmax < abs )
                {
                absmax  = abs ;
                imax    = i ;
                }
            }

        int signref = 0<cpref[imax] ? 1 : -1 ;

        //  the k'th hyperplane
        SEXP        svec = VECTOR_ELT(shyperplane,k) ;
        const int   *vec = INTEGER( svec ) ;
        int         nvec = Rf_length( svec );   //  nvec is the number of points in the hyperplane

        for( int j1=0 ; j1<nvec-1 ; j1++ )
            {
            int i = idxfromgnd[ vec[j1] ] ;     // 1-based

            for( int j2=j1+1 ; j2<nvec ; j2++ )
                {
                int j = idxfromgnd[ vec[j2] ] ; // 1-based

                int colidx  = PAIRINDEX( i, j, n ) ;    //  1-based

                double  *cp = crossprods + 3*(colidx-1) ; //  1-based to 0-based

                int sign =  0<cp[imax] ? 1 : -1 ;

                if( sign == signref )
                    {
                    cp[0]   = cpref[0] ;
                    cp[1]   = cpref[1] ;
                    cp[2]   = cpref[2] ;
                    }
                else
                    {
                    cp[0]   = -cpref[0] ;
                    cp[1]   = -cpref[1] ;
                    cp[2]   = -cpref[2] ;
                    }
                }
            }
        }

    R_Free( idxfromgnd );

    SEXP    out = PROTECT( Rf_allocVector(LGLSXP,1) );

    *(LOGICAL(out)) = true ;

    UNPROTECT(1);   //  variable out

    return(out);
    }



//  sgenidx     integer vector of 1-based indexes into a 3xn matrix of generators - smatgen.
//              These generators correspond to points in the hyperplane of a matroid.
//              These generators are coplanar in plane P.  They span a zonogon facet
//              and form the edges of this zonogon Z.  Usually there are only 2 of these.
//              The first generator in sgenidx is the edge from which the diameter is taken.
//              The diameter is from the center of this edge to the center of the antipodal edge.
//
//  snormal     a normal vector to the plane P, which gives an orientation to P and Z.
//
//  smatgen     3xn matrix of generators
//
//  scrossprods n(n-1)/2 matrix as returned from allcrossproducts(), possibly with normalization

SEXP
diametervector( SEXP sgenidx, SEXP snormal, SEXP smatgen, SEXP scrossprods )
    {
    int numgen  = Rf_length(sgenidx);

    if( numgen <= 1 )    return(R_NilValue);

    if( Rf_length(snormal) != 3 )    return(R_NilValue);

    const   int *dim    = INTEGER(Rf_getAttrib(smatgen,R_DimSymbol));

    int nrow    = dim[0] ;
    int n       = dim[1] ;

    if( nrow != 3 ) return(R_NilValue);

    dim = INTEGER(Rf_getAttrib(scrossprods, R_DimSymbol));
    if( dim[0] != 3 ) return(R_NilValue);
    if( dim[1] != (n*(n-1))/2 ) return(R_NilValue);

    //  find the largest coordinate of snormal
    int     kmax = -1;
    const   double  *normal = REAL(snormal);
    double  absmax=0 ;

    for( int k=0 ; k<3 ; k++ )
        {
        double  abs = fabs(normal[k]);
        if( absmax < abs )
            {
            absmax  = abs ;
            kmax    = k ;
            }
        }

    if( absmax == 0 )   return(R_NilValue);

    SEXP    out = PROTECT( Rf_allocVector(REALSXP,3) );
    double  *vecout = REAL(out);
    vecout[0]=vecout[1]=vecout[2]=0;

    const   int     *genidx     = INTEGER(sgenidx);
    const   double  *matgen     = REAL(smatgen);
    const   double  *crossprods = REAL(scrossprods);

    //  since numgen is usually 2, this for loop is usually executed only once.  optimize later ?

    for( int i=1 ; i<numgen ; i++ )
        {
#if 1
        if( genidx[i]<1  ||  n<genidx[i] )
            {
            Rprintf( "Internal Error.  genidx[%d]=%d.\n",  i, genidx[i] );
            UNPROTECT(1);   // variable out
            return(R_NilValue);
            }
#endif

        //  compare normal to the crossproduct
        int     diffidx = genidx[0] - genidx[i] ;
        int     pairidx;

        if( 0 < diffidx )
            pairidx =   PAIRINDEX( genidx[i], genidx[0], n ) - 1;   // 1-based to 0-based
        else
            pairidx =   PAIRINDEX( genidx[0], genidx[i], n ) - 1;   // 1-based to 0-based

#if 1
        if( pairidx<0  ||  (n*(n-1))/2 <= pairidx )
            {
            Rprintf( "Internal Error.  pairidx=%d.  genidx[0]=%d  genidx[%d]=%d.\n",
                        pairidx, genidx[0], i, genidx[i] );
            UNPROTECT(1);   // variable out
            return(R_NilValue);
            }
#endif

        double  sign;
        sign    = diffidx * normal[kmax] * crossprods[kmax + pairidx*3] ;
        sign    = 0<sign  ?  1  :  -1 ;

        for( int k=0 ; k<3 ; k++ )
            vecout[k] += sign * matgen[k + (genidx[i]-1)*3] ;   // 1-based to 0-based
        }

    UNPROTECT(1);   // variable out

    return( out );
    }



//  vec     numeric vector
//  n       length of vec, n>=2 not checked
//  skip    integer coordinate to skip (ignore)

bool
allequalskip( const double vec[], int n, int skip )
    {
    double  value   = vec[(skip+1) % n] ;

    for( int k=0 ; k<n ; k++ )
        {
        if( k == skip ) continue ;

        if( vec[k] != value )   return(false);
        }

    return(true);
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
//      *)  3 x m matrix with vector diameters in the columns
//      *)  logical m-vector with i'th value TRUE meaning that the point 0 is in the i'th facet

SEXP
diametermatrix( SEXP shyper, SEXP shypersub, SEXP scube, SEXP sgen, SEXP sground, SEXP snormal, SEXP smatgen, SEXP scrossprods )
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

    double  *pcube = R_Calloc( n, double );
    memcpy( pcube, REAL(scube), n*sizeof(*pcube) );




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

    SEXP    smat = PROTECT( Rf_allocMatrix(REALSXP,3,m) );
    double  *matout = REAL(smat);
    memset( matout, 0, 3*m*sizeof(*matout) ) ;

    SEXP    sfacet0 = PROTECT( Rf_allocVector(LGLSXP,m) ) ;
    int     *facet0 = LOGICAL(sfacet0) ;
    memset( facet0, 0, m*sizeof(*facet0) ) ;


    for( int i=0 ; i<m ; i++ )
        {
        //  get the i'th hyperplane.
        int         hyperidx = hypersub[i] ;

        if( hyperidx<1  ||  numhypers<hyperidx )
            {
            Rprintf( "Internal Error. hyperidx=%d  is invalid.  numhypers=%d.\n",
                        hyperidx, numhypers );
            UNPROTECT(2);   // variables smat, sfacet0
            return(R_NilValue);
            }

        SEXP        svec = VECTOR_ELT(shyper,hyperidx-1) ;   //  convert from 1-based to 0-based
        const int   *vec = INTEGER( svec ) ;
        int         nvec = Rf_length( svec );   //  nvec is the number of points in the hyperplane

        //  get the i'th normal
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
            UNPROTECT(2);   // variables smat, sfacet0
            return(R_NilValue);
            }

        double  signmultiplier  = (i==0) ? 0.5 : 1.0 ;

        //  get pointer to the i'th column of the output matrix
        double  *vecout = matout + i*3 ;

        bool    match = false;

        for( int k=0 ; k<nvec ; k++ )
            {
            //  get the k'th point in the hyperplane and convert to raw 1-based index genkidx
            int     genkidx = idxfromgnd[ vec[k] ];

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
                continue ;
                }

#if 1
            if( pairidx<0  ||  (n*(n-1))/2 <= pairidx )
                {
                Rprintf( "Internal Error. i=%d  k=%d  pairidx=%d.  genidx=%d  genkidx=%d.\n",
                            i, k, pairidx, genidx, genkidx );
                UNPROTECT(2);   // variables smat, sfacet0
                return(R_NilValue);
                }
#endif


            double  sign;
            sign    = diffidx * normal[jmax] * crossprods[jmax + pairidx*3] ;
            sign    = 0<sign  ?  1  :  -1 ;

            for( int j=0 ; j<3 ; j++ )
                vecout[j] += sign * matgen[j + (genkidx-1)*3] ;   // 1-based to 0-based

            //  advance to next center of edge
            pcube[genkidx-1] += signmultiplier * sign ;
            }

        if( ! match )
            {
            Rprintf( "Internal Error. i=%d  distinguished point %d not found in hyperplane %d.\n",
                        i, gen, hyperidx );
            UNPROTECT(2);   // variables smat, sfacet0
            return(R_NilValue);
            }


#if 0
        //  verify that all pcube values are +-0.5 except genidx, which is 0
        int badcount=0 ;
        for( int j=0 ; j<n ; j++ )
            {
            bool    ok = (j==genidx-1) ? (pcube[j]==0) : (fabs(pcube[j])==0.5) ;
            if( ! ok )  badcount++ ;
            }

        if( 0 < badcount )
            {
            Rprintf( "Internal Error. i=%d  pcube has %d bad values.\n", i, badcount );
            UNPROTECT(2);   // variables smat, sfacet0
            return(R_NilValue);
            }
#endif

        if( allequalskip( pcube, n, genidx-1 ) )
            {
            facet0[  i  ]       = true ;    //  add current facet to facet0
            facet0[ (i+1)%m ]   = true ;    //  add *next* facet to facet0 too !  Might be redundant and remove this later.
            }
        }


    R_Free( pcube );
    R_Free( idxfromgnd );

    SEXP    out = PROTECT( Rf_allocVector(VECSXP,2) );
    SET_VECTOR_ELT( out, 0, smat );
    SET_VECTOR_ELT( out, 1, sfacet0 );

    UNPROTECT(3);   // variables smat, sfacet0, and  out

    return( out );
    }





int
veccount( SEXP x )
    {
    switch( TYPEOF(x) )
        {
        case INTSXP :
            return 1 ;
            break ;

        case VECSXP :
            return Rf_length(x) ;
            break ;
        }

    return 0 ;
    }

bool
fillbuff( SEXP x, SEXP buff[], int *k )
    {
    int m ;

    switch( TYPEOF(x) )
        {
        case INTSXP :
            buff[ (*k)++ ] = x ;
            break;

        case VECSXP :
            m = Rf_length(x);
            for( int i=0 ; i<m ; i++ )  buff[ (*k)++ ] = VECTOR_ELT(x,i);
            break;
        }

    return true ;
    }



int
maxover( SEXP x )
    {
    int     *vec, n;

    int     out = 0;

    switch( TYPEOF(x) )
        {
        case INTSXP :
            vec = INTEGER(x);
            n   = Rf_length(x);
            for( int k=0 ; k<n ; k++ )
                {
                //Rprintf( "loop[%d]=%d    gmax=%d\n", k, loop[k], gmax );

                if( vec[k] < 1 )
                    {
                    Rprintf( "maxover(). ERR.  loop[%d] = %d < 1, which is invalid.\n", k, vec[k] );
                    return( 0 );
                    }

                out    = MAX2( vec[k], out ) ;
                }
            break;

        case VECSXP :
            n   = Rf_length(x);
            for( int k=0 ; k<n ; k++ )
                {
                SEXP        svec = VECTOR_ELT(x,k) ;
                const int   *vec = INTEGER( svec ) ;
                int         nvec = Rf_length( svec );

                for( int i=0 ; i<nvec ; i++ )
                    {
                    if( vec[i] < 1 )
                        {
                        Rprintf( "maxover(). ERR.  vec[%d] = %d < 1, which is invalid.\n", i, vec[i] );
                        return( 0 );
                        }

                    out    = MAX2( vec[i], out ) ;
                    }
                }
            break ;
        }

    return( out );
    }

//  xN  one of:
//      *) a list of integer vectors
//      *) an integer vector
//      *) NULL
//  all integer vectors are interpreted as sets of small positive integers
//
//      returns the union of all inputs, in strictly ascending order
SEXP
fastunion( SEXP x1, SEXP x2, SEXP x3 )
    {

    //Rprintf( "typeof x1 = %d\n", TYPEOF(x1) );
    //Rprintf( "typeof x2 = %d\n", TYPEOF(x2) );
    //Rprintf( "typeof x3 = %d\n", TYPEOF(x3) );


    int n = 0 ;
    n += veccount( x1 );
    n += veccount( x2 );
    n += veccount( x3 );

    //Rprintf( "allocating for %d vectors.  size=%d\n", n, sizeof(SEXP) );

    if( n == 0 )    return R_NilValue ;

    SEXP    *buff = R_Calloc( n, SEXP );

    int k = 0;

    fillbuff( x1, buff, &k );
    fillbuff( x2, buff, &k );
    fillbuff( x3, buff, &k );

    if( k != n )
        {
        Rprintf( "fastunion(). ERR. internal error %d != %d\n", k, n  );
        return( R_NilValue );
        }

    //  Pass 1.  find the maximum integer in all the sets
    int m;
    int imax    = 0;
    for( k=0 ; k<n ; k++ )
        {
        int *vec = INTEGER( buff[k] );
        m   = Rf_length( buff[k] );

        for( int j=0 ; j<m ; j++ )  imax = MAX2( vec[j], imax );
        }

    //Rprintf( "imax=%d\n", imax );

    int     *mask = R_Calloc( imax+1, int );

    //  Pass 2.  mark each position in mask[] with a 1
    for( k=0 ; k<n ; k++ )
        {
        int *vec = INTEGER( buff[k] );
        m   = Rf_length( buff[k] );

        for( int j=0 ; j<m ; j++ )  mask[ vec[j] ] = 1;
        }

    R_Free( buff );

    //  count the number of 1s
    int count = 0;
    for( int i=0 ; i<=imax ; i++ )   count += mask[i] ;

    SEXP    out = PROTECT( Rf_allocVector(INTSXP,count) );

    int     *pout   = INTEGER(out) ;

    k   = 0;
    for( int i=0 ; i<=imax ; i++ )
        {
        if( mask[i] )   pout[k++] = i ;
        }

    R_Free( mask );

    UNPROTECT(1);

    return( out );
    }


//  shyper      a list of integer vectors, each one of them defining a subset, and in increasing order.  NOT verified.
//  sground     an integer vector - the union of all the sets in shyper - and in increasing order.  NOT verified
//  sloop       an integer vector, which is disjoint from sground, this is verified.
//  smultiple   a list of integer vectors, each of them intersects sground in a single point, this is verified.
//              all the integers in smultiple must be disjoint from sloop; this is NOT verified.
//
//  returns a list the same length as shyper, with each integer vector augmented as appropriate by sloop and smultiple
//
//  all integers in the vectors are positive

SEXP
unsimplify( SEXP shyper, SEXP sground, SEXP sloop, SEXP smultiple )
    {
    const int     *ground = INTEGER(sground);
    int     nground = Rf_length(sground);
    int     gmax = ground[ nground-1 ];

    const int     *loop = INTEGER(sloop);
    int     nloop = Rf_length(sloop);

    int     nmultiple = Rf_length(smultiple);

    //  iterate over sloop and smultiple to find new gmax
    //  and test for negative values at the same time

    int     max1 = maxover( sloop );
    int     max2 = maxover( smultiple );

    gmax    = MAX2( max1, gmax );
    gmax    = MAX2( max2, gmax );


    //Rprintf( "gmax=%d\n", gmax );

    //  make a mask for the ground set
    int     *setmask = R_Calloc( gmax+1, int );

    for( int k=0 ; k<nground  ; k++ )
        setmask[ ground[k] ] = 1;


    //  check that sloop is disjoint from sground

    for( int k=0 ; k<nloop ; k++ )
        {
        if( setmask[ loop[k] ] )
            {
            Rprintf( "unsimplify(). ERR.  Point %d is in both ground and loop.\n", loop[k] );
            R_Free( setmask );
            return( R_NilValue );
            }
        }

    //  for each set in multiple, compute intersection with ground
    //  SEXP     *multiple = VECTOR_PTR(smultiple);

    int     *inter = R_Calloc( nmultiple, int );

    for( int k=0 ; k<nmultiple ; k++ )
        {
        SEXP        svec = VECTOR_ELT(smultiple,k) ;
        const int   *vec = INTEGER( svec ) ;
        int         nvec = Rf_length( svec );

        int     count = 0 ;

        for( int i=0 ; i<nvec ; i++ )
            {
            if( setmask[ vec[i] ] == 0 )    continue ;

            if( ++count == 1 )
                {
                inter[k] = i;
                //Rprintf( "inter[%d] = %d\n", k, inter[k] );
                }
            else
                {
                Rprintf( "unsimplify(). ERR.  Intersection of multiple #%d and ground set is not a singleton.\n", k+1 );
                R_Free( setmask );
                R_Free( inter );
                return( R_NilValue );
                }
            }

        if( count == 0 )
            {
            Rprintf( "unsimplify(). ERR.  Intersection of multiple %d and ground set is empty.\n", k+1 );
            R_Free( setmask );
            R_Free( inter );
            return( R_NilValue );
            }
        }


    int     nhyper = Rf_length(shyper);

    SEXP    out = PROTECT( Rf_allocVector(VECSXP,nhyper) );

    for( int i=0 ; i<nhyper ; i++ )
        {
        memset( setmask, 0, (gmax+1)*sizeof(int) );

        int     count=0 ;   //  number in the output set
        int     imax = 0;   //  max integer in the output set

        //  start with the i'th hyperplane itself
        SEXP        svec = VECTOR_ELT(shyper,i) ;
        const int   *vec = INTEGER( svec ) ;
        int         nvec = Rf_length( svec );

        for( int k=0 ; k<nvec ; k++ )
            {
            setmask[ vec[k] ] = 1;
            imax    = MAX2(vec[k],imax);
            }
        count   += nvec ;

        //  add the loops
        for( int k=0 ; k<nloop ; k++ )
            {
            setmask[ loop[k] ] = 1 ;
            imax    = MAX2(loop[k],imax);
            }
        count   += nloop ;


        //  add the extra multiples
        for( int k=0 ; k<nmultiple ; k++ )
            {
            SEXP        svec = VECTOR_ELT(smultiple,k) ;
            const int   *vec = INTEGER( svec ) ;

            if( setmask[ vec[ inter[k] ] ]  )
                {
                //  add the *other* multiples to setmask[]
                int         nvec = Rf_length( svec );
                for( int j=0 ; j<nvec ; j++ )
                    {
                    setmask[ vec[j] ] = 1 ;
                    imax    = MAX2(vec[j],imax);
                    }

                count   += nvec-1 ;     //  -1 because vec[ inter[k] ]  has already been marked 1
                }
            }

        if( count == Rf_length( svec ) )
            {
            //  a special case
            //  nothing was added to svec, so all we have to do is duplicate it
            //  NB:  this implies that there are no loops
            SET_VECTOR_ELT( out, i, Rf_duplicate(svec) );
            //  SET_VECTOR_ELT( out, i, svec );     this seems to work OK, but seems dangerous
            continue ;
            }

        SEXP    hpout = PROTECT( Rf_allocVector(INTSXP,count) );
        int     *vecout = INTEGER(hpout);

        int     m=0 ;
        for( int j=1 ; j<=imax ; j++ )
            {
            if( setmask[j] )    vecout[m++] = j ;
            }

        if( m != count )
            {
            Rprintf( "unsimplify().  ERR.  Internal %d != %d.\n", m, count );
            R_Free( inter );
            R_Free( setmask );
            UNPROTECT(2) ;      //  for hpout and out
            return( R_NilValue );
            }

        SET_VECTOR_ELT( out, i, hpout );

        UNPROTECT(1);   //  for the variable hpout
        }

    R_Free( inter );
    R_Free( setmask );

    UNPROTECT(1);    // for the variable out

    return( out );
    }




//  shyper      a list of integer vectors, each one of them defining a subset, and in increasing order.  NOT verified.
//              if not in increasing order, the output list may not be increasing too.
//  sground     an integer vector - the union of all the sets in shyper - and in increasing order.  NOT verified
//  sloop       an integer vector, which is a subset of sground, this is NOT verified.
//  smultiple   a list of integer vectors, each of them is a subset of sground,  NOT verified
//              all the integers in smultiple must be disjoint from sloop; this is NOT verified.
//
//  returns a list the same length as shyper, with sloop and the all but the first of smultiple removed
//
//  all integers in the vectors are positive
//
//  since shyper is a list of hyperplanes, if one point of a multiple group is present, then they all are,
//  and in particular the smallest (first) point.  This means we can remove all points, except the first one.
//  Compare with simplifygeneral().

SEXP
simplify( SEXP shyper, SEXP sground, SEXP sloop, SEXP smultiple )
    {
    const int     *ground = INTEGER(sground);
    int     nground = Rf_length(sground);
    int     gmax = ground[ nground-1 ];

    Rbyte   *maskremove = R_Calloc( gmax+1, Rbyte );

    const int     *loop = INTEGER(sloop);
    int     nloop = Rf_length(sloop);

    //  mark all loops for removal
    for( int k=0 ; k<nloop ; k++ )  { maskremove[ loop[k] ] = 1 ; }

    int     nmultiple = Rf_length(smultiple);

    for( int k=0 ; k<nmultiple ; k++ )
        {
        SEXP        svec = VECTOR_ELT(smultiple,k) ;
        const int   *vec = INTEGER( svec ) ;
        int         nvec = Rf_length( svec );

        //  in the next loop, note that vec[0] is skipped
        //  so this first or smallest point will NOT be removed
        for( int i=1 ; i<nvec ; i++ )   { maskremove[ vec[i] ] = 1 ; }
        }

    int     nhyper = Rf_length(shyper);

    SEXP    out = PROTECT( Rf_allocVector(VECSXP,nhyper) );

    for( int k=0 ; k<nhyper ; k++ )
        {
        //  get the k'th hyperplane
        SEXP        svec = VECTOR_ELT(shyper,k) ;
        const int   *vec = INTEGER( svec ) ;
        int         nvec = Rf_length( svec );

        //  Pass #1. count the number of points in output set
        int     count=0 ;   //  number in the output set; always a subset of vec[]
        for( int i=0 ; i<nvec ; i++ )
            {
            if( ! maskremove[ vec[i] ] ) { count++ ; }
            }

        if( count == nvec )
            {
            //  a special case
            //  nothing was removed from svec, so all we have to do is duplicate it
            //  NB:  this implies that there are no loops
            SET_VECTOR_ELT( out, k, Rf_duplicate(svec) );
            //  SET_VECTOR_ELT( out, k, svec );     this seems to work OK, but seems dangerous
            continue ;
            }

        //  Pass #2. load the output set, using exactly the same test
        SEXP    hpout = PROTECT( Rf_allocVector(INTSXP,count) );
        int     *vecout = INTEGER(hpout);

        int     m=0 ;
        for( int i=0 ; i<nvec ; i++ )
            {
            if( ! maskremove[ vec[i] ] ) { vecout[m++] = vec[i] ; }
            }

        SET_VECTOR_ELT( out, k, hpout );

        UNPROTECT(1);   //  for the variable hpout
        }

    R_Free( maskremove );

    UNPROTECT( 1 );    // for the variable out

    return( out );
    }



//  slist       a list of integer vectors, each one of them defining a subset of sground, and in increasing order.  NOT verified.
//              if not in increasing order, the output list may not be increasing too.
//  sground     an integer vector - the union of all the sets in slist - and in increasing order.  NOT verified
//  sloop       an integer vector, which is a subset of sground, this is NOT verified.
//  smultiple   a list of integer vectors, each of them is a subset of sground,  NOT verified
//              all the integers in smultiple must be disjoint from sloop; this is NOT verified.
//
//  returns a list the same length as slist, with sloop and the all but the first of smultiple removed
//
//
//  all integers in the vectors are positive
//
//  slist can have general subsets of ground, and they are NOT assumed to be hyperplanes.
//  All points of a specific multiple group are replaced by the first point, even if that point is not present.
//  simplifygeneral() uses more scratch vectors [5] than simplify() [1].
//  If a vector in slist has a point that is out of range, then a 0-length vector is returned.

SEXP
simplifygeneral( SEXP slist, SEXP sground, SEXP sloop, SEXP smultiple )
    {
    const int     *ground = INTEGER(sground);
    int     nground = Rf_length(sground);
    int     gmax = ground[ nground-1 ];


    //  maskloop marks all loops.  It does not change.
    Rbyte   *maskloop = R_Calloc( gmax+1, Rbyte );
    const int     *loop = INTEGER(sloop);
    int     nloop = Rf_length(sloop);

    //  mark all loops for removal
    for( int k=0 ; k<nloop ; k++ )  { maskloop[ loop[k] ] = 1 ; }

    //  idxgroup records the 1-based index of the multiple groups.
    //  test for no duplicates at the same time.
    int     *idxgroup = R_Calloc( gmax+1, int );

    int     nmultiple = Rf_length(smultiple);

    //  mingroup records the minimum point of each group
    int     *mingroup = R_Calloc( nmultiple+1, int );     // add +1 to make the group index 1-based

    for( int k=0 ; k<nmultiple ; k++ )
        {
        SEXP        svec = VECTOR_ELT(smultiple,k) ;
        const int   *vec = INTEGER( svec ) ;
        int         nvec = Rf_length( svec );

        mingroup[k+1] = INT_MAX ;

        for( int i=0 ; i<nvec ; i++ )
            {
            int     p = vec[i];

            if( maskloop[ p ] )
                {
                //  group k contains a loop, which is forbidden
                Rprintf( "simplifygeneral(). Internal error.  group %d contains a loop %d.\n", k, p );
                R_Free( maskloop );
                R_Free( idxgroup );
                R_Free( mingroup );
                return( R_NilValue );
                }

            if( idxgroup[ p ] )
                {
                //  group k intersects another group, which is forbidden
                Rprintf( "simplifygeneral(). Internal error.  group %d intersects group %d.  point %d\n",
                                    k, idxgroup[p]-1, p );
                R_Free( maskloop );
                R_Free( idxgroup );
                R_Free( mingroup );
                return( R_NilValue );
                }

            //  record the index of the group here
            idxgroup[p]  = k+1 ;     //  add 1 to make it a 1-based index

            mingroup[k+1]   = MIN2( p, mingroup[k+1] );
            }
        }

    //  firstgroup[] records whether the group has been encountered when scanning each set in slist
    bool    *metgroup = R_Calloc( nmultiple+1, bool );

    //  vecset[] is a scratch buffer to hold the modified set from slist
    int     *vecset = R_Calloc( gmax, int );     // add +1 to make the group index 1-based

    int     nlist = Rf_length(slist);

    SEXP    out = PROTECT( Rf_allocVector(VECSXP,nlist) );

    for( int k=0 ; k<nlist ; k++ )
        {
        //  get the k'th set
        SEXP        svec = VECTOR_ELT(slist,k) ;
        const int   *vec = INTEGER( svec ) ;
        int         nvec = Rf_length( svec );

        memset( metgroup, 0, (nmultiple+1)*sizeof(*metgroup) );

        int     count=0 ;

        for( int i=0 ; i<nvec ; i++ )
            {
            int     p = vec[i];

            if( p<1  ||  gmax<p )   { count=0 ; break ; }   //  maybe set count=-1 ?

            if( maskloop[p] )   continue ;  //  p is a loop so ignore it

            int idx = idxgroup[p] ;     //  1-based

            if( idx == 0 )
                {
                //  p is not in a group, so just add it
                vecset[count++] = p ;
                }
            else if( ! metgroup[idx] )
                {
                //  this group has not been encountered yet,
                //  so output the group minimum,
                //  and mark the group as being met, so no other points in the group are output
                vecset[count++] = mingroup[idx] ;
                metgroup[idx]   = true ;
                }
            }

        SEXP    setout = PROTECT( Rf_allocVector(INTSXP,count) );
        int     *vecout = INTEGER(setout);

        for( int i=0 ; i<count ; i++ )
            vecout[i] = vecset[i] ;

        SET_VECTOR_ELT( out, k, setout );

        UNPROTECT(1);   //  for the variable setout
        }

    R_Free( vecset );
    R_Free( maskloop );
    R_Free( idxgroup );
    R_Free( mingroup );
    R_Free( metgroup );

    UNPROTECT( 1 );    // for the variable out

    return( out );
    }


//  spair   Mx2 integer matrix of (i,j) pairs - with i < j
//              locating elements in the upper triangular part of an NxN matrix. these are 1-based
//  sn      the integer N
//
//  returns an integer M-vector; for each pair (i,j) returns the 1-based index of the pair
//          if i or j is NA_INTEGER, then the returned index is NA_INTEGER

SEXP
pairindex( SEXP spair, SEXP sn )
    {
    const   int *dim    = INTEGER(Rf_getAttrib(spair,R_DimSymbol)) ;

    int m   = dim[0] ;
    if( dim[1] != 2 )   return(R_NilValue);

    const   int *pairmat    = INTEGER(spair);

    int n   = *INTEGER(sn);


    SEXP    out = PROTECT( Rf_allocVector(INTSXP,m) );
    int     *idxvec = INTEGER(out);

    for( int k=0 ; k<m ; k++ )
        {
        //  int idx;
        int i = pairmat[k] ;        //  column #1
        int j = pairmat[k + m] ;    //  column #2

        idxvec[k]   = NA_INTEGER;

        if( i==NA_INTEGER  ||  j==NA_INTEGER )  continue ;

        if( 1<=i  &&  i<j  &&  j<=n )
            idxvec[k] = PAIRINDEX(i,j,n);
        }

    UNPROTECT(1);       //  variable out

    return( out );
    }






//  shyper      a list of integer vectors, each one of them defining a nontrivial hyperplane subset, and in increasing order.  Not checked
//  sground     an integer vector - the union of all the sets in shyper - and in increasing order.  Not checked
//
//
//  returns a list of 2-point hyperplanes that, combined with shyper, form a 2-partition of the ground set
//          if this is not possible, returns a list with info on the problem
//
//  Note:  The function really has 2 purposes:
//          1.  complete the given hyperplanes to a full set that satisfies the matroid axioms
//          2.  check that the given hyperplanes *can* be completed to a full set. Otherwise it's an ERROR.
//
//  all integers in the vectors are positive

SEXP
trivialhypers2( SEXP shyper, SEXP sground )
    {
    const int     *ground = INTEGER(sground);
    int     nground = Rf_length(sground);
    int     gmax = ground[ nground-1 ];

    //  make lookup table from ground points to raw index, 1..nground
    int     *idxfromgnd = R_Calloc( gmax+1, int );

    for( int k=0 ; k<nground ; k++ )
        idxfromgnd[ ground[k] ] = k+1 ;     //  this *was* just k

    //  make buffer to count the pairs
    int     npairs = (nground*(nground-1)) / 2 ;
    Rbyte   *paircount = R_Calloc( npairs, Rbyte );

    int     nhyper = Rf_length(shyper);

    Rbyte   cmax=1;
    int     csum=0;             //  the number of pairs contained in the given non-trivial hyperplanes
    int     pmax[2] = {-1,-1};  //  initialized to avoid a warning; anything will do

    for( int k=0 ; k<nhyper ; k++ )
        {
        SEXP    svec = VECTOR_ELT(shyper,k) ;

        const int   *vec = INTEGER( svec ) ;

        int     m = Rf_length(svec) ;

        if( m < 2 ) continue ;      // silently ignore

        for( int i=0 ; i<m-1 ; i++ )
            {
            int     ix   = idxfromgnd[ vec[i] ] ;   //  ix is 1-based
            //int     kdx0 = (jx*(jx-1)) / 2 ;

            for( int j=i+1 ; j<m ; j++ )
                {
                int jx = idxfromgnd[ vec[j] ];      //  jx is 1-based

                int kdx = PAIRINDEX(ix,jx,nground) - 1;   //  1-based to 0-based

                paircount[kdx] += 1;

                if( cmax < paircount[kdx] )
                    {
                    cmax = paircount[kdx] ;
                    pmax[0] = vec[i] ;
                    pmax[1] = vec[j] ;
                    }
                }
            }

        csum += (m*(m-1)) / 2 ;
        }

    if( 1 < cmax )
        {
        //  ! ERROR !
        R_Free( paircount );
        R_Free( idxfromgnd );

        //  return enough details for a good error message from the caller

        //  Rprintf( "trivialhypers2().  ERR. paircount[%d,%d] = %d.\n", pmax[0], pmax[1], cmax );
        SEXP    out = PROTECT( Rf_allocVector(VECSXP,2) );

        SEXP    svec ;

        svec = PROTECT( Rf_allocVector(INTSXP,1) );
        INTEGER(svec)[0] = cmax ;
        SET_VECTOR_ELT( out, 0, svec );

        svec = PROTECT( Rf_allocVector(INTSXP,2) );
        INTEGER(svec)[0] = pmax[0] ;
        INTEGER(svec)[1] = pmax[1] ;
        SET_VECTOR_ELT( out, 1, svec );

        UNPROTECT(2);   //  svec x 2

        SEXP    name2 = PROTECT(Rf_allocVector(STRSXP,2));
        SET_STRING_ELT(name2, 0, Rf_mkChar("cmax") );
        SET_STRING_ELT(name2, 1, Rf_mkChar("pmax") );
        Rf_setAttrib( out, R_NamesSymbol, name2 );
        UNPROTECT(1) ;  // name2

        UNPROTECT(1);   //  out

        return( out );
        }

    int     outcount    = npairs - csum ;

    if( outcount < 0 )
        {
        //  ERROR
        Rprintf( "trivialhypers2().  Internal Error. outcount = %d.\n", outcount );
        R_Free( paircount );
        R_Free( idxfromgnd );
        return( R_NilValue );
        }

    SEXP    out = PROTECT( Rf_allocVector(VECSXP,outcount) );

    int row=0 ;

    for( int i=1 ; i<nground ; i++  )
        {
        //  int     k0 = (j*(j-1)) / 2 ;

        for( int j=i+1 ; j<=nground ; j++  )
            {
            int k   = PAIRINDEX(i,j,nground) - 1;   //  1-based to 0-based

            if( paircount[k] == 0 )
                {
                SEXP    svec = PROTECT( Rf_allocVector(INTSXP,2) );
                int     *vec = INTEGER( svec );

                vec[0]  = ground[i-1];      //  back to 0-based
                vec[1]  = ground[j-1];      //  back to 0-based

                SET_VECTOR_ELT( out, row, svec );

                UNPROTECT(1);   // the variable svec
                //pout[ row ]             = ground[i] ;   //  column #1
                //pout[ row + outcount ]  = ground[j] ;   //  column #2
                row++ ;
                }
            }
        }


    R_Free( paircount );
    R_Free( idxfromgnd );

    UNPROTECT(1);    // the variable out

    if( row != outcount )
        {
        Rprintf( "trivialhypers2().  ERR.  Internal %d != %d.\n", row, outcount );
        return( R_NilValue );
        }

    return( out );
    }


//  sx      matrix of integers,  prevalidated
//  smargin an integer, 1 or 2, prevalidated

//  returns a list with items the rows (margin=1) or the columns (margin=2) of x

SEXP
matrix2list( SEXP sx, SEXP smargin )
    {
    matdat  md = extractmatdat( sx, smargin );
    if( md.imat == NULL )
        {
        return(R_NilValue);
        }

    SEXP  out = PROTECT( Rf_allocVector(VECSXP,md.nVec) );

    for( int j=0 ; j<md.nVec ; j++ )
        {
        int     *vec = md.imat + j*md.vecStep ;

        SEXP    svec = PROTECT( Rf_allocVector(INTSXP,md.vecLen) );

        int     *ovec = INTEGER(svec);

        for( int i=0 ; i<md.vecLen ; i++ )
            {
            ovec[i] = vec[i*md.eltStep] ;
            }

        SET_VECTOR_ELT( out, j, svec );

        UNPROTECT(1);   // variable svec
        }

    UNPROTECT(1);   //  variable out

    return(out);
    }


//  ssetlist    a list of integer vectors with no dups
//  sset        integer vector with no dups, representing a set

//  returns a logical vector out with out[i]=TRUE iff the i'th set is a subset of sset.
//  so length(out) = length(ssetlist)

SEXP
issubset( SEXP ssetlist, SEXP sset )
    {
    const int   *set = INTEGER(sset);
    int         nset = Rf_length(sset);

    //  pass #1 to find the max in set
    int     setmax = 0;
    for( int k=0 ; k<nset ; k++ )   { setmax = MAX2( set[k], setmax ) ; }

    //  pass #2 allocate mask and assign it
    Rbyte   *setmask = R_Calloc( setmax+1, Rbyte );
    for( int k=0 ; k<nset ; k++ )   { setmask[ set[k] ] = 1 ; }

    int     m = Rf_length(ssetlist) ;

    SEXP    out = PROTECT( Rf_allocVector(LGLSXP,m) );
    int     *pout = LOGICAL(out);

    for( int j=0 ; j<m ; j++ )
        {
        SEXP    svec = VECTOR_ELT(ssetlist,j);
        int     n = Rf_length(svec);

        if( nset < n )
            //  the j'th set is too big to be a subset of sset
            continue;

        const int   *setj = INTEGER( svec );
        int   i=0;
        for( i=0 ; i<n ; i++ )
            {
            if( setj[i]<1  ||  setmax<setj[i] ) break ;          //  setj[i] is too small or too big to be in sset

            if( setmask[ setj[i] ] == 0 )   break ; //  setj[i] is NOT in sset
            }

        pout[j] = (i == n) ;
        }

    R_Free( setmask ) ;

    UNPROTECT(1) ;      // variable out

    return(out);
    }



//  ssetlist    a list of integer vectors with no dups
//  sset        integer vector with no dups, representing a set

//  returns a logical vector out with out[i]=TRUE iff the i'th set contains sset.
//  so length(out) = length(ssetlist)

SEXP
issuperset( SEXP ssetlist, SEXP sset )
    {
    const int   *set = INTEGER(sset);
    int         nset = Rf_length(sset);

    //  pass #1 to find the max in set
    int     setmax = 0;
    for( int k=0 ; k<nset ; k++ )   { setmax = MAX2( set[k], setmax ) ; }

    //  pass #2 allocate mask and assign it
    Rbyte   *setmask = R_Calloc( setmax+1, Rbyte );
    for( int k=0 ; k<nset ; k++ )   { setmask[ set[k] ] = 1 ; }

    int     m = Rf_length(ssetlist) ;

    //  allocate output, but no need to initialize it
    SEXP    out = PROTECT( Rf_allocVector(LGLSXP,m) );
    int     *pout = LOGICAL(out);

    for( int j=0 ; j<m ; j++ )
        {
        SEXP    svec = VECTOR_ELT(ssetlist,j);
        int     n = Rf_length(svec);

        if( n < nset )
            {
            //  this set is too small to contain sset
            pout[j] = false ;
            continue;
            }

        //  count the number of points in set that setj contains
        int     count = 0 ;

        const int   *setj = INTEGER( svec );
        for( int i=0 ; i<n ; i++ )
            {
            if( setmax < setj[i] ) continue ;      //  setj[i] is too big

            count   += setmask[ setj[i] ] ;
            }

        pout[j] = (count == nset) ;
        }

    R_Free( setmask ) ;

    UNPROTECT(1) ;      // variable out

    return(out);
    }




//  ssetlist    a list of M integer vectors, with no dups, representing M sets
//  sset        integer vector with no dups, representing a set
//  sdecreasing a logical where TRUE means that the M sets are decreasing in size
//              when TRUE the function can be optimized, when FALSE it has no effect

//  returns a single logical value out with out=TRUE iff one of the M sets contains sset.

SEXP
anyissuperset( SEXP ssetlist, SEXP sset, SEXP sdecreasing )
    {
    const int   *set = INTEGER(sset);
    int         nset = Rf_length(sset);

    //  pass #1 to find the max in set
    int     setmax = 0;
    for( int k=0 ; k<nset ; k++ )   { setmax = MAX2( set[k], setmax ) ; }

    //  pass #2 allocate mask and assign it
    Rbyte   *setmask = R_Calloc( setmax+1, Rbyte );
    for( int k=0 ; k<nset ; k++ )   { setmask[ set[k] ] = 1 ; }

    const bool  decreasing = *(LOGICAL(sdecreasing)) ;

    int     m = Rf_length(ssetlist) ;

    //  allocate output
    SEXP    out = PROTECT( Rf_allocVector(LGLSXP,1) );
    int     *pout = LOGICAL(out);
    *pout   = false ;

    for( int j=0 ; j<m ; j++ )
        {
        SEXP    svec = VECTOR_ELT(ssetlist,j);
        int     n = Rf_length(svec);

        if( n < nset )
            {
            //  the j'th set is too small to contain sset
            if( decreasing )
                //  subsequent sets are too small as well, so quit
                break;
            else
                //  keep searching
                continue;
            }

        //  count the number of points in set that setj contains
        int     count = 0 ;

        const int   *setj = INTEGER( svec );
        for( int i=0 ; i<n ; i++ )
            {
            if( setmax < setj[i] ) continue ;      //  setj[i] is too big

            count   += setmask[ setj[i] ] ;

            if( count == nset ) break ;     //  count cannot get any larger, so quit
            }

        if( count == nset )
            {
            //  the j'th set contains the given set, so can assign *pout and quit
            *pout   = true;
            break ;
            }
        }

    R_Free( setmask ) ;

    UNPROTECT(1) ;      // variable out

    return(out);
    }







//  shyper      a list of integer vectors, each one of them defining a nontrivial hyperplane subset, and in increasing order.  Not checked
//  sground     an integer vector - the union of all the sets in shyper - and in increasing order.  Not checked
//
//  all integers in the vectors are positive
//
//  returns a list of two items
//      *) incident[] an integer vector.  incident[i] = # of hyperplanes that contain point i
//      *) hash[] a floating-point vector of point i, hash[i] can be very large, so use double

SEXP
incidencedata( SEXP shyper, SEXP sground )
    {
    const int     *ground = INTEGER(sground);
    int     nground = Rf_length(sground);
    int     gmax = ground[ nground-1 ];

    SEXP    out = PROTECT( Rf_allocVector(VECSXP,2) );

    SEXP    sincid      = PROTECT( Rf_allocVector(INTSXP,gmax) );
    int     *incident   = INTEGER(sincid);
    memset( incident, 0, gmax*sizeof(*incident) );

    SEXP    shash   = PROTECT( Rf_allocVector(REALSXP,gmax) );
    double  *hash   = REAL(shash);
    memset( hash, 0, gmax*sizeof(*hash) );

    SET_VECTOR_ELT( out, 0, sincid );
    SET_VECTOR_ELT( out, 1, shash );

    int     m = Rf_length(shyper) ;

    for( int j=0 ; j<m ; j++ )
        {
        SEXP    svec = VECTOR_ELT(shyper,j);
        const int   *hp = INTEGER( svec );
        int     n = Rf_length(svec);

        for( int k=0 ; k<n ; k++ )
            {
            incident[ hp[k]-1 ]++ ;     //  1-based to 0-based

            //  use 1.0 here to prevent signed integer overflow,
            //  which was detected as an Additional Issue May 29, 2023
            hash[ hp[k]-1 ] += (j+1.0)*(j+1.0)*M_LN2;      //  1-based to 0-based
            }
        }

    UNPROTECT(2);   //  variables sincid, shash

    //  name both items in the output list
    SEXP    name2 = PROTECT(Rf_allocVector(STRSXP,2));
    SET_STRING_ELT(name2, 0, Rf_mkChar("incident") );
    SET_STRING_ELT(name2, 1, Rf_mkChar("hash") );
    Rf_setAttrib( out, R_NamesSymbol, name2 );

    UNPROTECT(2);   //  variables out, name2

    return(out);
    }

//  shyper      a list of integer vectors, each one of them defining a nontrivial hyperplane subset, and in increasing order.  Not checked
//  sground     an integer vector - the union of all the sets in shyper - and in increasing order.  Not checked
//  ssubset     an integer vector, defining a subset of ground
//
//  all integers in the vectors are positive
//
//  returns a logical matrix that is length(shyper) x length(subset)
//  out[i,j] = TRUE  iff  the j'th element of ssubset is in the i'th hyperplane

SEXP
incidencematrix( SEXP shyper, SEXP sground, SEXP ssubset )
    {
    const int     *ground = INTEGER(sground);
    int     nground = Rf_length(sground);
    int     gmax = ground[ nground-1 ];

    //  make small lookup table, with 1-based index of ssubset
    int     *index  = R_Calloc( gmax+1, int );

    const int   *subset = INTEGER(ssubset);
    int     n = Rf_length(ssubset) ;
    for( int j=0 ; j<n ; j++ )
        index[ subset[j] ]  = j+1 ;     //  1-based index in subset, add 1 to distinguish from 0 and subtract later

    int     m = Rf_length(shyper) ;

    SEXP    out = PROTECT( Rf_allocMatrix(LGLSXP,m,n) );
    int     *matout = LOGICAL(out);

    memset( matout, 0, m*n*sizeof(*matout) ) ;

    for( int i=0 ; i<m ; i++ )
        {
        SEXP        svec = VECTOR_ELT( shyper, i );
        const int   *vec = INTEGER(svec);
        int         nhper = Rf_length(svec);

        for( int k=0 ; k<nhper ; k++ )
            {
            int     jdx = index[ vec[k] ];
            if( 0 < jdx )
                matout[ i + (jdx-1)*m ] = 1;   //  convert from 1-based to 0-based.  now subtracting 1
            }
        }

    R_Free( index );

    UNPROTECT(1);   //  out

    return(out);
    }



//  shyper      a list of integer vectors, each one of them defining a nontrivial hyperplane subset, and in increasing order.  Not checked
//              it also must be for a simple rank 3 matroid, which guarantees that each point is in at most n-1 hyperplanes.  Not checked
//  sground     an integer vector of length n - the union of all the sets in shyper - and in increasing order.  Not checked

//  all integers in the vectors are positive
//
//  returns an integer matrix that is n x (n-1)
//  The i'th row of the matrix contains the indexes of the hyperplanes that contain the i'th point in sground = the i'th belt of the zonohedron
//  The indexes are in increasing order, and do NOT match the order on the zonohedron.
//  If there are fewer than n-1 hyperplanes, the row is 0-terminated.
//  i is 1-based, i.e. in 1..n

SEXP
beltmatrix( SEXP shyper, SEXP sground )
    {
    int     nhyper = Rf_length(shyper);

    const int     *ground = INTEGER(sground);
    int     n   = Rf_length(sground);
    int     gmax = ground[ n-1 ];

    //  make lookup table from ground points to raw index, 0..n-1
    int     *idxfromgnd = R_Calloc( gmax+1, int );

    for( int k=0 ; k<n ; k++ )
        idxfromgnd[ ground[k] ] = k ;

    SEXP    out = PROTECT( Rf_allocMatrix(INTSXP,n,n-1) );
    int     *matout = INTEGER(out);
    memset( matout, 0, Rf_length(out)*sizeof(*matout) ) ;


    //  make vector to count the number of hyperplanes containing point i, so far
    int *count  = R_Calloc( n, int );

    bool    ok = true ;

    for( int k=0 ; k<nhyper && ok ; k++ )
        {
        SEXP    svec = VECTOR_ELT(shyper,k) ;
        const int   *vec = INTEGER( svec ) ;
        int         nvec = Rf_length( svec );

        for( int j=0 ; j<nvec ; j++ )
            {
            int     p = vec[j] ;            //  get the point in the ground set
            int     idx = idxfromgnd[p] ;   //  convert from ground set to 0-based index

            if( count[idx] == n-1 )
                {
                Rprintf( "beltmatrix().  Internal Error. count[%d] = %d.", idx, n-1 );
                ok  = false ;
                break;
                }

#if 0
            if( Rf_length(out) <= idx + count[idx]*n )
                {
                Rprintf( "beltmatrix().  Internal Error. idx=%d  count[idx]=%d.", idx, count[idx] );
                ok  = false ;
                break;
                }
#endif

            matout[ idx + count[idx]*n ] = k+1 ;    // k is 0-based, so convert to 1-based
            count[idx]++ ;
            }
        }

    R_Free( count );
    R_Free( idxfromgnd );

    UNPROTECT(1);   //  out

    if( ! ok )    return(R_NilValue);

    return(out);
    }


SEXP
duplicateR( SEXP x )
    {
    return Rf_duplicate( x );
    }


//  x       any R variable
//
//  returns either 1 or 2 4-byte integers equivalent to the address of x.
//  1 or 2 depending on whether it is compiled to 32-bit or 64-bit.

SEXP
obj_addr( SEXP x )
    {
    //  struct pair  { union { SEXP s; int i[2]; } ; };     the -Wpedantic option does *NOT* like this
    //  matroid.c:1910:48: warning: ISO C99 doesn't support unnamed structs/unions [-Wpedantic]

#if SIZEOF_SIZE_T == 8
    typedef union { SEXP s; int i[2]; } pair ;      //  the -Wpedantic option *does* like this

    pair p;

    p.s  = x;

    SEXP    out = PROTECT( Rf_allocVector(INTSXP,2) );

    INTEGER(out)[0] = p.i[0] ;
    INTEGER(out)[1] = p.i[1] ;

#elif SIZEOF_SIZE_T == 4
    typedef union { SEXP s; int i[1]; } pair ;      //  the -Wpedantic option *does* like this

    pair p;

    p.s  = x;

    SEXP    out = PROTECT( Rf_allocVector(INTSXP,1) );

    INTEGER(out)[0] = p.i[0] ;

#else
    #error "SIZEOF_SIZE_T is neither 8 nor 4."
#endif

    UNPROTECT(1);   //  out

    return out;
    }