
#include "R_ext/Print.h"

//  #include "config.h"         // this file is created by configure.win or configure


#define HAVE_CXX1X      // avoid warning  "Specified C++11: please drop specification unless essential"

#ifdef HAVE_CXX1X

#include <climits>  /* for CHAR_BIT */
#include <cstdint>   /* for uint_fast32_t etc */
#include <cstddef>  /* for size_t */

#define RANDBIT8  0x8B
#define RANDBIT16 0x135E
#define RANDBIT32 0x183B2DEA
#define RANDBIT64 0x3A0A3F58E697236C
union {
    uint_fast8_t ui8;
    uint_fast16_t ui16;
    uint_fast32_t ui32;
    uint_fast64_t ui64;
    size_t randbit;
} randbit;
unsigned int lshift, rshift;
extern "C" int initHash(void)
{
    switch (sizeof(size_t) * CHAR_BIT){
        case 32:    randbit.ui32 = RANDBIT32; lshift=6;  rshift=2; return 1;
        case 64:    randbit.ui64 = RANDBIT64; lshift=12; rshift=4; return 1;
        case 16:    randbit.ui16 = RANDBIT16; lshift=3;  rshift=1; return 1;
        case 8:     randbit.ui8  = RANDBIT8;  lshift=2;  rshift=0; return 1;
        default:    randbit.ui64 = RANDBIT64; lshift=12; rshift=4; return 0;
    }
    Rprintf("initHash() !\n" );
}
#include "rcSetHash.h"



// instantiation of global objects:
//  this only works because R is single-threaded  [GD]
vecSetHash<int>             intVecSetHash;
vecSetHash<double>          doubleVecSetHash;
vecSetHash<CharSEXP>        charsexpVecSetHash;
vecSetHash<Rcomplex>        cmplxVecSetHash;
vecSetHash<unsigned char>   rawVecSetHash;      // Rbyte is an alias of unsigned char


//  grpDuplicatedMat() uses these objects, and only works because R is single-threaded  [GD]
vecMapHash<int>             intVecMapHash;
vecMapHash<double>          doubleVecMapHash;
vecMapHash<CharSEXP>        charsexpVecMapHash;
vecMapHash<Rcomplex>        cmplxVecMapHash;
vecMapHash<unsigned char>   rawVecMapHash;      // Rbyte is an alias of unsigned char

extern "C" {

SEXP dupAtomMatHash(SEXP x, SEXP MARGIN, SEXP fromLast)
{/* returns a logical vector of duplicated rows of numeric matrix x */
    SEXP out;
    int* dim;
    dim=INTEGER(Rf_getAttrib(x, R_DimSymbol));
    out = PROTECT(Rf_allocVector(LGLSXP, dim[*INTEGER(MARGIN)-1]));

    switch (TYPEOF(x)) {
        case REALSXP:
            doubleVecSetHash.duplicatedMat  (REAL(x), dim, dim+1,  LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
            break;
        case INTSXP:  // factor type is also covered here
            // if(!inherits(x, "factor"))
                intVecSetHash.duplicatedMat (INTEGER(x), dim, dim+1,  LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
            // else {;}
            break;
        case LGLSXP:
            intVecSetHash.duplicatedMat (LOGICAL(x), dim, dim+1,  LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
            break;
        case STRSXP: {
            CharSEXP* charSexpPtr = new CharSEXP [ dim[0]*dim[1] ];
            for(int i=dim[0]*dim[1]-1; i>=0; --i)
                charSexpPtr[i].sexp = STRING_ELT(x, i);

            charsexpVecSetHash.duplicatedMat    (charSexpPtr, dim, dim+1, LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );

            delete[] charSexpPtr;
            break;
        }
        case CPLXSXP:
            cmplxVecSetHash.duplicatedMat   (COMPLEX(x), dim, dim+1,  LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
            break;
        case RAWSXP:
            rawVecSetHash.duplicatedMat (RAW(x), dim, dim+1,  LOGICAL(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
            break;
        default:
            Rf_error("C function 'dupAtomMatHash' only accepts REALSXP, LGLSXP, INTSXP and STRSXP");
    }

    UNPROTECT(1);
    return out;
}



SEXP anyDupAtomMatHash(SEXP x, SEXP MARGIN, SEXP fromLast)
{/* returns a logical vector of duplicated rows of numeric matrix x */
    SEXP out;
    int* dim;
    dim=INTEGER(Rf_getAttrib(x, R_DimSymbol));
    out = PROTECT(Rf_allocVector(INTSXP, 1));

    switch (TYPEOF(x)) {
        case REALSXP:
            doubleVecSetHash.anyDuplicatedMat   (REAL(x), dim, dim+1,  INTEGER(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
            break;
        case INTSXP:  // factor type is also covered here
            // if(!inherits(x, "factor"))
                intVecSetHash.anyDuplicatedMat  (INTEGER(x), dim, dim+1,  INTEGER(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
            // else {;}
            break;
        case LGLSXP:
            intVecSetHash.anyDuplicatedMat  (LOGICAL(x), dim, dim+1,  INTEGER(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
            break;
        case STRSXP: {
            CharSEXP* charSexpPtr = new CharSEXP [ dim[0]*dim[1] ];
            for(int i=dim[0]*dim[1]-1; i>=0; --i)
                charSexpPtr[i].sexp = STRING_ELT(x, i);

            charsexpVecSetHash.anyDuplicatedMat (charSexpPtr, dim, dim+1, INTEGER(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );

            delete[] charSexpPtr;
            break;
        }
        case CPLXSXP:
            cmplxVecSetHash.anyDuplicatedMat    (COMPLEX(x), dim, dim+1,  INTEGER(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
            break;
        case RAWSXP:
            rawVecSetHash.anyDuplicatedMat  (RAW(x), dim, dim+1,  INTEGER(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
            break;
        default:
            Rf_error("C function 'anyDupAtomMatHash' only accepts REALSXP, LGLSXP, INTSXP and STRSXP");
    }

    UNPROTECT(1);
    return out;
}


/* returns an integer vector of duplicated rows of numeric matrix x */
SEXP grpDupAtomMatHash( SEXP x, SEXP MARGIN )
    {
    int     *dim = INTEGER(Rf_getAttrib(x, R_DimSymbol));
    SEXP    out = PROTECT( Rf_allocVector(INTSXP, dim[*INTEGER(MARGIN)-1]) );

    bool    ok;
    int     nGrps[3];

    switch (TYPEOF(x)) {
        case REALSXP:
            ok = doubleVecMapHash.grpDuplicatedMat(REAL(x), dim, dim+1, *INTEGER(MARGIN)==1,  INTEGER(out), nGrps  );
            break;
        case INTSXP:  // factor type is also covered here
            // if(!inherits(x, "factor"))
                ok = intVecMapHash.grpDuplicatedMat(INTEGER(x), dim, dim+1, *INTEGER(MARGIN)==1,  INTEGER(out), nGrps  );
            // else {;}
            break;
        case LGLSXP:
            ok = intVecMapHash.grpDuplicatedMat(LOGICAL(x), dim, dim+1, *INTEGER(MARGIN)==1,  INTEGER(out), nGrps );
            break;
        case STRSXP: {
            CharSEXP* charSexpPtr = new CharSEXP [ dim[0]*dim[1] ];
            for(int i=dim[0]*dim[1]-1; i>=0; --i)
                charSexpPtr[i].sexp = STRING_ELT(x, i);

            ok = charsexpVecMapHash.grpDuplicatedMat(charSexpPtr, dim, dim+1, *INTEGER(MARGIN)==1, INTEGER(out), nGrps );

            delete[] charSexpPtr;
            break;
        }
        case CPLXSXP:
            ok = cmplxVecMapHash.grpDuplicatedMat(COMPLEX(x), dim, dim+1, *INTEGER(MARGIN)==1,  INTEGER(out), nGrps );
            break;
        case RAWSXP:
            ok = rawVecMapHash.grpDuplicatedMat(RAW(x), dim, dim+1, *INTEGER(MARGIN)==1,  INTEGER(out), nGrps );
            break;
        default:
            Rf_error( "C function 'grpDupAtomMatHash' only accepts REALSXP, LGLSXP, INTSXP, STRSXP, CPLXSXP, and RAWSXP.");
    }

    if( ! ok )
        Rf_error( "C function 'grpDupAtomMatHash' unexpected failure !" );

    SEXP nGroups;
    nGroups = PROTECT( Rf_allocVector(INTSXP,3) );
    INTEGER(nGroups)[0] = nGrps[0];
    INTEGER(nGroups)[1] = nGrps[1];
    INTEGER(nGroups)[2] = nGrps[2];

    Rf_setAttrib(out, Rf_install("ngroups"), nGroups);

    UNPROTECT(2);   //  out and nGroups

    return out;
    }

}

#endif
