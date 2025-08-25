
#include <R.h>
#include <Rinternals.h>

#include "macros.h"

//  #define TEST_RCOND


SEXP
trans2vertex( SEXP sn )
    {
    int     n = *INTEGER(sn);

    if( n <= 0 )
        {
        return(R_NilValue);
        }

    int rowsout = n*(n-1)+2 ;

    SEXP    out = PROTECT( Rf_allocMatrix(INTSXP,rowsout,2) );
    int     *matout = INTEGER(out);

    //  south pole
    matout[ 0 + 0*rowsout ] = 0;
    matout[ 0 + 1*rowsout ] = NA_INTEGER;

    int k   = 1 ;
    for( int count=1 ; count<n ; count++ )
        {
        for( int start=1 ; start<=n ; start++, k++ )
            {
            matout[ k + 0*rowsout ] = count;
            matout[ k + 1*rowsout ] = start;
            }
        }

    //  north pole
    matout[ k + 0*rowsout ] = n;
    matout[ k + 1*rowsout ] = NA_INTEGER;

    UNPROTECT(1);   //  variable out

    return(out);
    }



//  parameters:
//  N       the dimension of the cube, must be a positive integer
//  count   the number of 1s in the vertex, must be between 0 and N
//  start   the starting index of the run of 1s, 1-based
//
//  returns the row index in vertex[], as returned by trans2vertex()

static  int
vtxidxfromcode( int N, int count, int start )
    {
    if( count == 0 )    return( 1 ) ;

    if( count == N )    return( N*(N-1) + 2 ) ;

    return( N*(count-1) + start + 1 ) ;
    }



//  arguments:
//
//  N       dimension of the cube, must be a positive integer
//  scrange range for the count of 1s for the edges

SEXP
trans2edge( SEXP sn, SEXP scrange )
    {
    int     n = *INTEGER(sn);

    if( n <= 0 )    return(R_NilValue);

    if( n == 1 )
        {
        //  trivial case, only 1 edge, and only 2 vertices,  the 1-cube
        SEXP    out = PROTECT( Rf_allocMatrix(INTSXP,1,2) );
        INTEGER(out)[0] = 1 ;
        INTEGER(out)[1] = 2 ;
        UNPROTECT(1) ;   //  variable out
        return(out);
        }


    int *crange     = INTEGER(scrange);
    int countmin    = crange[0] ;
    int countmax    = crange[1];

    int midcountmin = MAX2( countmin, 1 );
    int midcountmax = MIN2( countmax, n-1 );

    //  count the
    int edges   = 0;
    if( countmin == 0 ) edges += n;

    if( midcountmin < midcountmax ) edges   += 2*n*(midcountmax - midcountmin) ;

    if( countmax == n ) edges += n;

    SEXP    out = PROTECT( Rf_allocMatrix(INTSXP,edges,2) );
    int     *matout = INTEGER(out);

    int i=0 ;       //  row index

    if( countmin == 0 )
        {
        //  add the "cap" at the south pole
        for( int start=1 ; start<=n ; start++, i++ )
            {
            matout[i + 0*edges] = vtxidxfromcode(n,0,NA_INTEGER);
            matout[i + 1*edges] = vtxidxfromcode(n,1,start) ;
            }
        }


    for( int count=midcountmin ; count<midcountmax ; count++ )
        {
        for( int start=1 ; start<=n ; start++ )
            {
            //  find index of current vertex
            int idx = vtxidxfromcode(n,count,start);

            //  add 1 on the left
            int sprev   = ((start-2 + n) % n) + 1 ;     //  start and sprev are 1-based

            matout[i + 0*edges] = idx ;
            matout[i + 1*edges] = vtxidxfromcode(n,count+1,sprev) ;
            i++ ;

            //  add 1 on the right
            matout[i + 0*edges] = idx ;
            matout[i + 1*edges] = vtxidxfromcode(n,count+1,start) ;
            i++ ;
            }
        }


    if( countmax == n )
        {
        //  add the "cap" at the north pole
        for( int start=1 ; start<=n ; start++, i++ )
            {
            matout[i + 0*edges] = vtxidxfromcode(n,n-1,start);
            matout[i + 1*edges] = vtxidxfromcode(n,n,NA_INTEGER) ;
            }
        }

    UNPROTECT(1) ;   //  variable out

    if( i != edges )
        {
        Rprintf( "trans2edge(). ERR. internal error %d != %d\n", i, edges  );
        return(R_NilValue);
        }

    return(out);
    }


//  sn      dimension of cube [-1/2,1/2]^N
//  scount  integer M-vector with the number of 1s in the vertex
//  sstart  integer M-vector with the starting index of the run of 1s, 1-based

//  returns an MxN matrix where the i'th row is the desired vertex of the N-cube, [-1/2,1/2]^N



SEXP
vertexfromcode( SEXP sn, SEXP scount, SEXP sstart )
    {
    int n   = *INTEGER(sn);

    int m   = Rf_length(scount);

    if( Rf_length(sstart) != m )    return(R_NilValue);

    const int *count  = INTEGER(scount);
    const int *start  = INTEGER(sstart);

    SEXP    out = PROTECT( Rf_allocMatrix(REALSXP,m,n) );
    double  *matout = REAL(out);

    //  initialize entire matrix to all -1/2
    for( int k=0 ; k<m*n ; k++ )
        matout[k] = -0.5 ;

    for( int i=0 ; i<m ; i++ )
        {
        if( count[i] == 0 ) continue ;      // south pole, all -1/2  already

        if( count[i] < n )
            {
            for( int k=0 ; k<count[i] ; k++ )
                {
                int j   = (start[i]-1 + k) % n ;    //  convert from 1-based to 0-based
                matout[i + j*m] = 0.5;
                }
            }
        else
            {
            //  count[i]==n  =>  north pole, all 1/2
            for( int j=0 ; j<n ; j++ )  { matout[i + j*m] = 0.5 ; }
            }
        }

    UNPROTECT(1) ;   //  variable out

    return( out );
    }


//  smatgen     M x N matrix of N generators, defining a 2-transition subcomplex / surface in R^M
//  sgensum     matrix of the same size - the cumsum over the columns of matgen
//
//  returns N*(N-1)/2 x M matrix the facet centers, in the same order as allpairs()
//          these facet centers are in *uncentered* surface coords

SEXP
allpgramcenters2trans( SEXP smatgen, SEXP sgensum )
    {
    const int *dim1  = INTEGER(Rf_getAttrib(smatgen, R_DimSymbol));
    int m   = dim1[0] ;
    int n   = dim1[1] ;

    const int *dim2  = INTEGER(Rf_getAttrib(sgensum, R_DimSymbol));

    if( dim2[0] != m  ||  dim2[1] != n )    return(R_NilValue);

    const   double  *matgen = REAL(smatgen);
    const   double  *gensum = REAL(sgensum);

    int     pgrams  = (n*(n-1)) / 2;

    SEXP    out = PROTECT( Rf_allocMatrix(REALSXP,pgrams,m) );
    double  *center = REAL(out);

    int     count   = 0;
    for( int j1=0 ; j1<n-1 ; j1++ )
        {
        const double    *v1 = matgen + m*j1;
        const double    *s1 = gensum + m*j1;

        for( int j2=j1+1 ; j2<n ; j2++, count++ )
            {
            const double    *v2 = matgen + m*j2;
            const double    *s2 = gensum + m*(j2-1);

            for( int k=0 ; k<m ; k++ )
                center[ count + k*pgrams ]  = 0.5 * ( v1[k] + v2[k] );

            if( j1 < j2-1 )
                {
                for( int k=0 ; k<m ; k++ )
                    center[ count + k*pgrams ] += s2[k] - s1[k] ;
                }
            }
        }

    UNPROTECT(1) ;   //  variable out

    return( out );
    }


//  scenterrot  N x 3 matrix with N pgram centers in the rows
//  sbaserot    3-vector of a basepoint of ray, with ray direction parallel to (0,0,1)  -  the z-axis
//
//  searches to find the row of the matrix whose first 2 coords are closest to those of sbaserot
//  however if the 3rd coord of the row is less than sbaserot[2], the row is ignored
//
//  returns the 0-based index of the optimal row

SEXP
optimalcenter( SEXP scenterrot, SEXP sbaserot )
    {
    const int *dim1  = INTEGER(Rf_getAttrib(scenterrot, R_DimSymbol));

    int n   = dim1[0] ;

    if( n == 0  ||  dim1[1] != 3 )  return(R_NilValue);

    if( Rf_length(sbaserot) != 3 )  return(R_NilValue);

    const   double  *centerrot = REAL(scenterrot);
    const   double  *baserot = REAL(sbaserot);

    int     imin = -1 ;
    double  dmin = FLT_MAX ;

    for( int i=0 ; i<n ; i++ )
        {
        double  z = centerrot[ i + 2*n ] ; // row=i col=3

        if( z <= baserot[2] )    continue ;  //  ignore this row

        double  x = centerrot[ i ]  -  baserot[0] ;
        double  y = centerrot[ i + n ] - baserot[1] ;

        double  d =  x*x + y*y ;
        if( d < dmin )
            {
            dmin    = d ;
            imin    = i ;
            }
        }

    if( imin < 0 )  return(R_NilValue) ;    //  no row found

    SEXP    out = PROTECT( Rf_allocVector(INTSXP,1) );

    *INTEGER(out)   = imin ;

    UNPROTECT(1) ;   //  variable out

#if 0
    double  x = centerrot[ imin ]  -  baserot[0] ;
    double  y = centerrot[ imin + n ] - baserot[1] ;
    double  z = centerrot[ imin + 2*n ] ; // row=i col=3
    Rprintf( "optimalcenter().  xyz=%g,%g,%g.  imin=%d\n", x, y, z, imin ) ;
#endif

    return( out );
    }




//  scenterrot  N x 3 matrix with N pgram centers in the rows
//  sbaserot    3-vector of a basepoint of ray, with ray direction parallel to (0,0,1)  -  the z-axis
//  sidxpair    N x 2 matrix of 1-based indexes, matching the rows of scenterrot
//  sgenrot     3 x M matrix of rotated generators
//
//  we must have N = M*(M-1)

//  returns a list with:
//      isect   the 0-based index of the row of scenterrot that contains the basepoint
//      alpha   2-vector of *centered* coefficents, both have abs value <= 1/2

SEXP
findpgram2D( SEXP scenterrot, SEXP sbaserot, SEXP sidxpair, SEXP sgenrot )
    {
    const int *dim ;

    dim = INTEGER(Rf_getAttrib(scenterrot, R_DimSymbol));

    int n   = dim[0] ;

    if( n == 0  ||  dim[1] != 3 )  return(R_NilValue);

    if( Rf_length(sbaserot) != 3 )  return(R_NilValue);

    dim = INTEGER(Rf_getAttrib(sidxpair, R_DimSymbol));

    if( dim[0] != n  ||  dim[1] != 2 )  return(R_NilValue);

    dim = INTEGER(Rf_getAttrib(sgenrot, R_DimSymbol));

    int m   = dim[1] ;

    if( dim[0] != 3  ||  m*(m-1) != n )  return(R_NilValue);


    const   double  *centerrot  = REAL(scenterrot);
    const   double  *baserot    = REAL(sbaserot);
    const   int     *idxpair    = INTEGER(sidxpair);
    const   double  *genrot     = REAL(sgenrot);

    int     isect   = -1 ;  //  index of pgram that intersects
    double  alphagood[2] = { NA_REAL, NA_REAL };

#ifdef  TEST_RCOND
    double  tol     = 1.e-14 ;
#endif

    for( int i=0 ; i<n ; i++ )
        {
        //  get 1-based indexes
        int     j1 = idxpair[i] ;
        int     j2 = idxpair[i + n] ;

#if 0
        int     diff    = abs(j1 - j2) ;    //  integer absolute value
        if( diff==1  ||  diff==(m-1) )
            //  these pgrams are already handled by the testpole() section
            continue ;
#endif

        const   double  *gen1 = genrot + 3*(j1-1) ;     //  1-based to 0-based
        const   double  *gen2 = genrot + 3*(j2-1) ;     //  1-based to 0-based

        //  get the z-coord of the center of the pgram
        double  z = centerrot[ i + 2*n ] ; // row=i col=3

        if( z  +  (fabs(gen1[2]) + fabs(gen2[2]))/2.0  <  baserot[2] )
            //  this pgram is below the plane of constant z, containing the basepoint
            //  so the ray cannot possibly intersect this pgram
            continue ;

        //  compute and test the determinant
        double  det = gen1[0] * gen2[1] - gen1[1] * gen2[0] ;
        if( det == 0 )  continue ;      //  degenerate pgram

#ifdef  TEST_RCOND
        //  compute and test reciprocal condition number, with Froebenius norm
        double  norm2 = SQUARE(gen1[0]) + SQUARE(gen1[1]) ;
        norm2   +=  SQUARE(gen2[0]) + SQUARE(gen2[1]) ;

        double  rcond   = fabs(det) / norm2 ;

        if( rcond <= tol )
            {
            Rprintf( "%s(). rcond=%g <= %g.  near-degenerate pgram %d skipped.\n", __func__, rcond, tol, i ) ;
            continue ;
            }
#endif

        double  x = baserot[0] - centerrot[ i ] ;
        double  y = baserot[1] - centerrot[ i + n ] ;

        //  use Cramer's rule to find alpha[]
        double  alpha[2] ;

        alpha[0]    = (x * gen2[1] - y * gen2[0] ) / det;

        if( 0.5 < fabs(alpha[0]) )  continue ;

        alpha[1]    = (gen1[0] * y - gen1[1] * x ) / det;

        //  Rprintf( "j1=%d  j2=%d  alpha = %g,%g\n", j1, j2, alpha[0], alpha[1] );

        //  if( 0.5 < fabs(alpha[0]) )  continue ;
        if( 0.5 < fabs(alpha[1]) )  continue ;

        //  compute z at intersection and test it
        double  zinter  = z + alpha[0]*gen1[2]  +  alpha[1]*gen2[2] ;

        if( zinter <= baserot[2] )   continue ;  //  the wrong side (negative) of the ray


        //  we have found it
        isect           = i;
        alphagood[0]    = alpha[0];
        alphagood[1]    = alpha[1];

        break ;
        }

    SEXP    out = PROTECT( Rf_allocVector(VECSXP,2) );

    SEXP    sisect      = PROTECT( Rf_allocVector(INTSXP,1) );
    *INTEGER(sisect)    = isect ;

    SEXP    salpha  = PROTECT( Rf_allocVector(REALSXP,2) );

    REAL(salpha)[0] = alphagood[0] ;
    REAL(salpha)[1] = alphagood[1] ;

    SET_VECTOR_ELT( out, 0, sisect );
    SET_VECTOR_ELT( out, 1, salpha );

    UNPROTECT(3);   //  sisect and salpha and variable out

    return( out );
    }

