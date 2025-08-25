
#include <stdint.h>         // for uint64_t  and  uint8_t
#include <R.h>
#include <Rinternals.h>
//  #include "R_ext/Print.h"

//  #include "config.h"         // this file is created by configure.win or configure

typedef uint64_t nanotime_t;


//  #define  DO_TIMES_EXECUTION

#ifdef  DO_TIMES_EXECUTION

#if defined(WIN32)
#include "nanotimer_windows.h"
#elif defined(MB_HAVE_MACH_TIME)
#include "nanotimer_macosx.h"
#elif defined(MB_HAVE_CLOCK_GETTIME) && defined(MB_CLOCKID_T)
#include "nanotimer_clock_gettime.h"
#elif defined(MB_HAVE_GETHRTIME)
#include "nanotimer_rtposix.h"
#elif defined(MB_HAVE_GETTIMEOFDAY)
#include "nanotimer_gettimeofday.h"
#else /* ./configure should prevent this, but just in case... */
#error "Unsupported OS."
#endif

#endif


//  #include <climits>  /* for CHAR_BIT */
//  #include <cstdint>   /* for uint_fast32_t etc */
//  #include <cstddef>  /* for size_t */

#include "flVecMap.h"



// instantiation of globals to force Rbyte to be compiled
typedef flVecMap<Rbyte>     rawflVecMap;      // Rbyte is an alias of unsigned char


//  the following struct assumes that ints are 32-bit and pointers are 64-bit
struct  pointer_converter
    {
private:
    typedef union { void *p; int i[2]; } ppair ;        //  the -Wpedantic option *does* like this

public:
    static  void *  pointerfromint( const int i[] )
                    {
                    ppair   thepair;

                    thepair.i[0] = i[0] ;
                    thepair.i[1] = i[1] ;

                    return( thepair.p );
                    }


    static  void    intfrompointer( const void *p, int i[] )
                    {
                    ppair   thepair;

                    thepair.p  = (void *) p;

                    i[0]    = thepair.i[0];
                    i[1]    = thepair.i[1];
                    }
    };



extern "C" {

//  n       length of all raw vectors
//  count   expected number of vectors in the map

//  makes a new rawflVecMap, and returns a pointer to it, as 2 integers

SEXP makeRawMap( SEXP sn, SEXP scount )
    {
    int     n       = *( INTEGER(sn) );
    int     count   = *( INTEGER(scount) );

    //  rawflVecMap     dummy(n,count) ;

    rawflVecMap *pmap    = new rawflVecMap(n,count) ;

    //Rprintf( "makeRawMap(). allocated pointer %p.\n", (void *)pmap );

    //  return pmap as 2 integers
    //ppair   thepair;
    //thepair.p  = pmap;

    SEXP    out = PROTECT( Rf_allocVector(INTSXP,2) );

    pointer_converter::intfrompointer( pmap, INTEGER(out) );

    //INTEGER(out)[0] = thepair.i[0] ;
    //INTEGER(out)[1] = thepair.i[1] ;

    UNPROTECT(1);       //  out

    return out;
    }


//  smap    two 4-byte integers, to be converted to a pointer to the rawflVecMap
//  svec    vector of raws, all values 0 or 1  (not verified).
//          length must be a multiple of pmap->getN()
//  scomp   boolean, take complement of svec
//
//  returns vector of indexes, with length equal to length(svec) / pmap->getN()

SEXP getIndexRaw( SEXP smap, SEXP svec, SEXP scomp )
    {
#ifdef  DO_TIMES_EXECUTION
    nanotime_t  time_start = get_nanotime();
#endif

    //  smap must be an integer vector with length 2
    bool    ok =  TYPEOF(smap)==INTSXP   &&   Rf_length(smap)==2 ;
    if( ! ok )
        {
        Rprintf( "getIndexRaw().  smap is invalid." );
        return( R_NilValue );
        }

    //  convert smap to a pointer
    //ppair   thepair;

    //thepair.i[0] = INTEGER(smap)[0] ;
    //thepair.i[1] = INTEGER(smap)[1] ;

    rawflVecMap *pmap   = (rawflVecMap *) pointer_converter::pointerfromint( INTEGER(smap) ) ;

    //Rprintf( "getIndexRaw(). recovered pointer %p.\n", (void *)pmap );

    if( pmap == NULL )
        {
        Rprintf( "getIndexRaw().  smap pointer is NULL !" );
        return( R_NilValue );
        }


    //  svec must be a raw vector with length a multiple of pmap->getN()
    int lenraw  = Rf_length(svec) ;

    int n       = int(pmap->getN()) ;

    ok =  TYPEOF(svec)==RAWSXP   &&  0<lenraw  &&  (lenraw % n) == 0;
    if( ! ok )
        {
        Rprintf( "getIndexRaw().  svec is invalid.  typeof=%d  length=%d.", TYPEOF(svec), Rf_length(svec) );
        return( R_NilValue );
        }

    //  scomp must be a single logical
    ok  = TYPEOF(scomp)==LGLSXP  &&  Rf_length(scomp)==1 ;
    if( ! ok )
        {
        Rprintf( "getIndexRaw().  scomp is invalid." );
        return( R_NilValue );
        }

    int     count = lenraw / n ;

    SEXP    out = PROTECT( Rf_allocVector(INTSXP,count) );

    for( int k=0 ; k<count ; k++ )
        {
        INTEGER(out)[k] = pmap->getIndex( RAW(svec) + k*n, n, LOGICAL(scomp)[0] ) ;
        }

    UNPROTECT(1);       //  out

#ifdef  DO_TIMES_EXECUTION
    pmap->mTimeInsertion  += get_nanotime() - time_start;
#endif

    return out;
    }


//  smap    two 4-byte integers, to be converted to a pointer to the rawflVecMap
//  smatgen 3xN matrix of doubles - the simplified generators of the zonohedron
//          N must be equal to pmap->getN()
//
//  returns 3xV matrix of doubles, where V is the number of binary codes in the map
//          and the number of vertices in the zonohedron
//
//  This function only uses additions, and avoids all multiplications.

SEXP computeVertices( SEXP smap, SEXP smatgen )
    {
#ifdef  DO_TIMES_EXECUTION
    nanotime_t  time_start = get_nanotime();
#endif

    //  smap must be an integer vector with length 2
    bool    ok ;

    ok  =  TYPEOF(smap)==INTSXP   &&   Rf_length(smap)==2 ;
    if( ! ok )
        {
        Rprintf( "computeVertices().  smap is invalid." );
        return( R_NilValue );
        }

    rawflVecMap *pmap   = (rawflVecMap *) pointer_converter::pointerfromint( INTEGER(smap) ) ;

    //Rprintf( "computeVertices(). recovered pointer %p.\n", (void *)pmap );

    if( pmap == NULL )
        {
        Rprintf( "computeVertices().  smap pointer is NULL !" );
        return( R_NilValue );
        }

    int n       = int(pmap->getN()) ;


    //  smatgen must be a 3xn matrix of doubles

    int *dim    = INTEGER(Rf_getAttrib(smatgen, R_DimSymbol));
    int nrow    = dim[0] ;
    int ncol    = dim[1] ;

    ok  =  TYPEOF(smatgen)==REALSXP   &&  nrow==3  &&  ncol==n ;
    if( ! ok )
        {
        Rprintf( "computeVertices().  smatgen is invalid." );
        return( R_NilValue );
        }

    const   double  *matgen = REAL(smatgen);


    const   std::map<  std::vector<uint8_t>, int> &themap = pmap->mMap ;

    int     verts = int( themap.size() );

    SEXP    out = PROTECT( Rf_allocMatrix(REALSXP,3,verts) );

    double  *outmat = REAL(out) ;

    std::map<  std::vector<uint8_t>, int >::const_iterator    p ;

    long double accum[3] ;

    for( p=themap.begin() ; p!=themap.end() ; p++ )
        {
        //  get reference to code[], whose values are 0 or 1, not verified

        const   std::vector<uint8_t>    &code   = p->first ;

        //  find the column in outmat[] to compute
        //  from the algorithm, this is the only time that this column is encountered
        int     col = p->second - 1 ;       //  1-based to 0-based column index

        if( verts <= col )
            {
            Rprintf( "computeVertices().  col is invalid.  col=%d  >=  %d=verts.", col, verts );
            UNPROTECT(1);   //  out
            return( R_NilValue );
            }

        //  compute column col in outmat[]
        int j = 3*col ;

        accum[0] = 0;
        accum[1] = 0;
        accum[2] = 0;

        for( int k=0 ; k<n ; k++ )
            {
            if( code[k] )
                {
                //  code[k] is not 0, and therefore must be 1
                //  increment the column by generator k
                accum[0] += matgen[3*k + 0] ;
                accum[1] += matgen[3*k + 1] ;
                accum[2] += matgen[3*k + 2] ;
                }
            }

        outmat[j+0] = accum[0];
        outmat[j+1] = accum[1];
        outmat[j+2] = accum[2];
        }

    UNPROTECT(1);   //  out

#ifdef  DO_TIMES_EXECUTION
    pmap->mTimeVertices = get_nanotime() - time_start;
#endif

    return( out );
    }


//  smap    two 4-byte integers, to be converted to a pointer to the rawflVecMap
//
//  returns NxV matrix of raw 0-1 values
//              N is the length of the raw codes in the map
//              V is the number of binary codes in the map
//          and the number of vertices in the zonohedron
//

SEXP getCodes( SEXP smap )
    {
    //  smap must be an integer vector with length 2
    bool    ok ;

    ok  =  TYPEOF(smap)==INTSXP   &&   Rf_length(smap)==2 ;
    if( ! ok )
        {
        Rprintf( "getCodes().  smap is invalid." );
        return( R_NilValue );
        }

    rawflVecMap *pmap   = (rawflVecMap *) pointer_converter::pointerfromint( INTEGER(smap) ) ;

    //Rprintf( "getCodes(). recovered pointer %p.\n", (void *)pmap );

    if( pmap == NULL )
        {
        Rprintf( "getCodes().  smap pointer is NULL !" );
        return( R_NilValue );
        }

    int n   = int(pmap->getN()) ;

    const   std::map< std::vector<uint8_t>, int > &themap = pmap->mMap ;

    int     verts = int( themap.size() );

    SEXP    out = PROTECT( Rf_allocMatrix(RAWSXP,n,verts) );

    uint8_t *outmat = RAW(out) ;

    ::memset( outmat, 255, n*verts ) ;  //  for debugging

    std::map<  std::vector<uint8_t>, int >::const_iterator    p ;

    for( p=themap.begin() ; p!=themap.end() ; p++ )
        {
        //  get reference to code[], whose values are 0 or 1, not verified

        const   std::vector<uint8_t>    &code   = p->first ;

        //  find the column in outmat[] to load
        //  from the algorithm, this is the only time that this column is encountered
        int     col = p->second - 1 ;       //  1-based to 0-based column index

        if( verts <= col )
            {
            Rprintf( "getCodes().  col is invalid.  col=%d  >=  %d=verts.", col, verts );
            UNPROTECT(1);   //  out
            return( R_NilValue );
            }

        ::memcpy( outmat + col*n, & code[0], n );
        }


    UNPROTECT(1);   //  out

    return( out );
    }



SEXP deleteRawMap( SEXP smap )
    {
    //  smap must be an integer vector of length 2
    bool    ok =  TYPEOF(smap)==INTSXP   &&   Rf_length(smap)==2 ;
    if( ! ok )
        {
        Rprintf( "getIndexRaw().  smap is invalid." );
        return( R_NilValue );
        }

    //  convert smap to a pointer
    //ppair   thepair;

    //thepair.i[0] = INTEGER(smap)[0] ;
    //thepair.i[1] = INTEGER(smap)[1] ;

    //rawflVecMap *pmap   = (rawflVecMap *) thepair.p ;

    rawflVecMap *pmap   = (rawflVecMap *) pointer_converter::pointerfromint( INTEGER(smap) ) ;

    //  Rprintf( "deleteRawMap(). recovered pointer %p.\n", (void *)pmap );

#ifdef  DO_TIMES_EXECUTION
    Rprintf( "deleteRawMap(). Total insertion time: %g sec.  Vertex computation time: %g sec.\n",
                    1.e-9 * pmap->mTimeInsertion, 1.e-9 * pmap->mTimeVertices );
#endif

    delete  pmap ;  pmap = NULL ;

    //  clear the input pointer, so the dangling pointer cannot cause any trouble
    INTEGER(smap)[0]    = 0;
    INTEGER(smap)[1]    = 0;

    SEXP    out = PROTECT( Rf_allocVector(LGLSXP,1) );

    LOGICAL(out)[0] = true ;

    UNPROTECT(1);       //  out

    return out;
    }



}       //  extern "C"
