

//  in clipquad3D() we clip a quadrangle against the positive octant
//  since there are 4 vertices input, and 3 planes,
//  there are at most 7 vertices output
    
//  vertout[] must be big enough to hold 7 vertices

#define VERTSOUTMAX 7

extern  bool    clipquad3D( double quad[3][4], double *vertout[3], int *nout );
