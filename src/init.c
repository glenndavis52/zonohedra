
//  #include "config.h"         // this file is created by configure.win or configure

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>



//	#define HASHED_NAME(fun) fun ## Hash
extern SEXP dupAtomMatHash(SEXP, SEXP, SEXP);
extern SEXP anyDupAtomMatHash(SEXP, SEXP, SEXP);
extern SEXP grpDupAtomMatHash(SEXP, SEXP );
extern int  initHash(void);

//	#define HASHED_NAME(fun) fun


extern  SEXP    dbl_dig(void);
extern  SEXP    collapseGroups1D_R( SEXP, SEXP ) ;
extern  SEXP    conditionalAntipodal( SEXP, SEXP, SEXP ) ;
extern  SEXP    normalizeMatrix( SEXP, SEXP ) ;
extern  SEXP    duplicateR( SEXP x ) ;
extern  SEXP    obj_addr( SEXP x ) ;
extern  SEXP    allcrossproducts( SEXP sx );
extern  SEXP    allpairs( SEXP sground );
extern  SEXP    fastunion( SEXP x1, SEXP x2, SEXP x3 ) ;
extern  SEXP    simplify( SEXP shyper, SEXP sground, SEXP sloop, SEXP smultiple );
extern  SEXP    simplifygeneral( SEXP slist, SEXP sground, SEXP sloop, SEXP smultiple );
extern  SEXP    unsimplify( SEXP shyper, SEXP sground, SEXP sloop, SEXP smultiple );
extern  SEXP    trivialhypers2( SEXP shyper, SEXP sground );
extern  SEXP    matrix2list( SEXP sx, SEXP smargin );
extern  SEXP    incidencedata( SEXP shyper, SEXP sground );
extern  SEXP    incidencematrix( SEXP shyper, SEXP sground, SEXP ssubset );
extern  SEXP    issubset( SEXP ssetlist, SEXP sset ) ;
extern  SEXP    issuperset( SEXP ssetlist, SEXP sset ) ;
extern  SEXP    anyissuperset( SEXP ssetlist, SEXP sset, SEXP sdecreasing );
extern  SEXP    whichMaxMatrix( SEXP sx, SEXP smargin ) ;
extern  SEXP    trans2vertex( SEXP sn ) ;
extern  SEXP    trans2edge( SEXP sn, SEXP scrange );
extern  SEXP    vertexfromcode( SEXP sn, SEXP scount, SEXP sstart );
extern  SEXP    pairindex( SEXP spair, SEXP sn );
extern  SEXP    beltmatrix( SEXP shyper, SEXP sground ) ;
extern  SEXP    cumsumMatrix( SEXP sx, SEXP smargin );
extern  SEXP    diametervector( SEXP sgenidx, SEXP snormal, SEXP smatgen, SEXP scrossprods );
extern  SEXP    diametermatrix( SEXP shyper, SEXP shypersub, SEXP scube, SEXP sgen, SEXP sground, SEXP snormal, SEXP smatgen, SEXP scrossprods );
extern  SEXP    multicopy( SEXP sdestmat, SEXP sdiff, SEXP ssrcmat, SEXP sdestidx );
extern  SEXP    beltdata( SEXP shyper, SEXP shypersub, SEXP scube, SEXP sgen, SEXP sground, SEXP snormal, SEXP smatgen, SEXP scrossprods );
extern  SEXP    radiusfacet( SEXP shyper, SEXP sground, SEXP sradiusgen );
extern  SEXP    sectionzonohedron( SEXP shyper, SEXP sfacetcenter, SEXP sfacetnormal, SEXP scenternormal, SEXP sbeta,
                        SEXP sground, SEXP sgennormal, SEXP smatgen, SEXP scrossprods ) ;
extern  SEXP    plusEqual( SEXP smat, SEXP svec, SEXP smargin );
extern  SEXP    timesEqual( SEXP smat, SEXP svec, SEXP smargin );
extern  SEXP    sumMatVec( SEXP smat, SEXP svec, SEXP smargin );
extern  SEXP    beltmidpoints( SEXP shyper, SEXP shypersub, SEXP sgen, SEXP scenter, SEXP snormal, SEXP sground, SEXP smatgen, SEXP scrossprods );
extern  SEXP    allpgramcenters2trans( SEXP smatgen, SEXP sgensum ) ;
extern  SEXP    snapcrossprods( SEXP scrossprods, SEXP shyperplane, SEXP scrossprodsref, SEXP sground );
extern  SEXP    area_sphtri( SEXP pa, SEXP pb, SEXP pc );
extern  SEXP    linkingnumber( SEXP smatgen, SEXP sidxpair, SEXP scenter, SEXP spoint );
extern  SEXP    linkingnumber2( SEXP smatcum, SEXP spoint );
extern  SEXP    linkingnumber3( SEXP smatgen, SEXP sidxpair, SEXP scenter, SEXP spoint );
extern  SEXP    optimalcenter( SEXP scenterrot, SEXP sbaserot );
extern  SEXP    extend_antipodal( SEXP smat );
extern  SEXP    findpgram2D( SEXP scenterrot, SEXP sbaserot, SEXP sidxpair, SEXP sgenrot );
extern  SEXP    transitioncount( SEXP spt );
extern  SEXP    rotation2pole_test( SEXP u );
extern  SEXP    dist2pgram_test( SEXP spoint, SEXP sv1, SEXP sv2, SEXP scenter, SEXP snormal );
extern  SEXP    dist2surface( SEXP smatgen, SEXP sidxpair, SEXP scenter, SEXP snormal, SEXP spoint );
extern  SEXP    clipquad( SEXP smatquad ) ;
extern  SEXP    makeRawMap( SEXP sn, SEXP scount );
extern  SEXP    getIndexRaw( SEXP smap, SEXP svec, SEXP scomp );
extern  SEXP    computeVertices( SEXP smap, SEXP smatgen );
extern  SEXP    getCodes( SEXP smap );
extern  SEXP    deleteRawMap( SEXP smap );

static R_CallMethodDef callMethods[]  = {
  {"dupAtomMatHash", (DL_FUNC) &dupAtomMatHash, 3},
  {"anyDupAtomMatHash", (DL_FUNC) &anyDupAtomMatHash, 3},
  {"grpDupAtomMatHash", (DL_FUNC) &grpDupAtomMatHash, 2},
  {"dbl_dig", (DL_FUNC) &dbl_dig, 0},
  {"collapseGroups1D_R", (DL_FUNC) &collapseGroups1D_R, 2},
  {"conditionalAntipodal", (DL_FUNC) &conditionalAntipodal, 3},
  {"normalizeMatrix", (DL_FUNC) &normalizeMatrix, 2},
  {"allcrossproducts", (DL_FUNC) &allcrossproducts, 1},
  {"allpairs", (DL_FUNC) &allpairs, 1},
  {"fastunion", (DL_FUNC) &fastunion, 3},
  {"simplify", (DL_FUNC) &simplify, 4},
  {"simplifygeneral", (DL_FUNC) &simplifygeneral, 4},
  {"unsimplify", (DL_FUNC) &unsimplify, 4},
  {"duplicateR", (DL_FUNC) &duplicateR, 1},
  {"trivialhypers2", (DL_FUNC) &trivialhypers2, 2},
  {"matrix2list", (DL_FUNC) &matrix2list, 2},
  {"incidencedata", (DL_FUNC) &incidencedata, 2},
  {"incidencematrix", (DL_FUNC) &incidencematrix, 3},
  {"issubset", (DL_FUNC) &issubset, 2},
  {"issuperset", (DL_FUNC) &issuperset, 2},
  {"anyissuperset", (DL_FUNC) &anyissuperset, 3},
  {"obj_addr", (DL_FUNC) &obj_addr, 1},
  {"whichMaxMatrix", (DL_FUNC) &whichMaxMatrix, 2},
  {"trans2vertex", (DL_FUNC) &trans2vertex, 1},
  {"trans2edge", (DL_FUNC) &trans2edge, 2},
  {"vertexfromcode", (DL_FUNC) &vertexfromcode, 3},
  {"pairindex", (DL_FUNC) &pairindex, 2},
  {"beltmatrix", (DL_FUNC) &beltmatrix, 2},
  {"cumsumMatrix", (DL_FUNC) &cumsumMatrix, 2},
  {"diametervector", (DL_FUNC) &diametervector, 4},
  {"diametermatrix", (DL_FUNC) &diametermatrix, 8},
  {"multicopy", (DL_FUNC) &multicopy, 4},
  {"beltdata", (DL_FUNC) &beltdata, 8},
  {"radiusfacet", (DL_FUNC) &radiusfacet, 3},
  {"sectionzonohedron", (DL_FUNC) &sectionzonohedron, 9},
  {"plusEqual", (DL_FUNC) &plusEqual, 3},
  {"timesEqual", (DL_FUNC) &timesEqual, 3},
  {"sumMatVec", (DL_FUNC) &sumMatVec, 3},
  {"beltmidpoints", (DL_FUNC) &beltmidpoints, 8},
  {"allpgramcenters2trans", (DL_FUNC) &allpgramcenters2trans, 2},
  {"snapcrossprods", (DL_FUNC) &snapcrossprods, 4},
  {"area_sphtri", (DL_FUNC) &area_sphtri, 3},
  {"linkingnumber", (DL_FUNC) &linkingnumber, 4},
  {"linkingnumber2", (DL_FUNC) &linkingnumber2, 2},
  {"linkingnumber3", (DL_FUNC) &linkingnumber3, 4},
  {"optimalcenter", (DL_FUNC) &optimalcenter, 2},
  {"extend_antipodal", (DL_FUNC) &extend_antipodal, 1},
  {"findpgram2D", (DL_FUNC) &findpgram2D, 4},
  {"transitioncount", (DL_FUNC) &transitioncount, 1},
  {"rotation2pole_test", (DL_FUNC) &rotation2pole_test, 1},
  {"dist2pgram_test", (DL_FUNC) &dist2pgram_test, 5},
  {"dist2surface", (DL_FUNC) &dist2surface, 5},
  {"clipquad", (DL_FUNC) &clipquad, 1},
  {"makeRawMap", (DL_FUNC) &makeRawMap, 2},
  {"getIndexRaw", (DL_FUNC) &getIndexRaw, 3},
  {"computeVertices", (DL_FUNC) &computeVertices, 2},
  {"getCodes", (DL_FUNC) &getCodes, 1},
  {"deleteRawMap", (DL_FUNC) &deleteRawMap, 1},
  {NULL, NULL, 0}
};


void R_init_zonohedra(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
   R_useDynamicSymbols(info, FALSE);
   R_forceSymbols(info, TRUE);

   if( !initHash() )    Rf_error("Hashing initialization error");
}
