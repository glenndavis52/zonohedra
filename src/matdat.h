

typedef struct
{
    double  *mat ;      // pointer to the matrix of doubles
    int     *imat ;     // pointer to the matrix of integers
    int     dim[2] ;    // rows and columns
    int     eltStep ;   // step between elements of x
    int     vecStep ;   // step between vectors of x
    int     vecLen ;    // length of each vector
    int     nVec ;      // number of vectors
} matdat;


extern  matdat  extractmatdat( SEXP sx, SEXP smargin ) ;

