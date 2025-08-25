
//      useful macros       //



#define MAX2(a, b) (((a) >= (b)) ? (a) : (b))
#define MIN2(a, b) (((a) <= (b)) ? (a) : (b))

#define ABS(a) (((0 <= (a)) ? (a) : -(a))

#define SQUARE(a)   ((a)*(a))

#define SIGNOF(x)   ( 0<(x) ? 1 : ( (x)<0 ? -1 : 0 ) )


//  i and j are row and column indexes, both 1-based.  The matrix is nxn.
//  The return value is also 1-based.
#define PAIRINDEX( i, j, n )    ( (i-1)*n - ((i)*(i+1))/2 + j )

