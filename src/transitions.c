
#include <stdbool.h>


#include <R.h>
#include <Rinternals.h>



SEXP
transitioncount( SEXP spt )
    {
    const   double  *pt = REAL(spt);

    int     n   = Rf_length(spt);

    //  make logical mask of the interior points
    //  and find the first boundary point
    bool    *interior = R_Calloc( n, bool );
    
    int     bfirst  = -1;
    
    for( int i=0 ; i<n ; i++ )
        {
        interior[i] = 0 < pt[i]  &&  pt[i] < 1 ;
        
        if( ! interior[i]  &&  bfirst<0 )   bfirst = i ;
        }
        

    if( bfirst < 0 )
        {
        //  all values are interior, a special case
        R_Free( interior );
        
        SEXP    out = PROTECT( Rf_allocVector(INTSXP,1) );
        *INTEGER(out)   =  2 * ( (n+1)/2 ) ;
        UNPROTECT(1) ;
        return( out ) ;
        }
        
    int     count = 0;
        
    //  initialize the first run
    double  valfirst = pt[bfirst] ;     //  0 or 1
    
    int     ifirst = bfirst ;
    
    for( int k=0 ; k<n ; k++  )
        {
        //  iprev   = (k   + bfirst) % n ;
        
        //  compute the index i into both pt[] and interior[]
        //  the first time through the loop, k=0 and i is one past bfirst
        //  the last time, k=n-1 and i=bfirst
        int     i       = (k+1 + bfirst) % n ;
        
        //  bool    inter_prev  = interior[iprev];
        //  bool    inter       = interior[i];
        
        if( interior[i] )
            count++ ;
        else
            {
            //  not an interior value, pt[i] is 0 or 1,  end of a run
            bool    equal   = valfirst == pt[i] ;           // are 'boundary endpoints' of the run are equal
            int     len     = (i-1 - ifirst + n) % n ;      // the number of interior values in the run
            if( ((int)equal) == (len % 2) )   count++ ;     // e.g. if len==0 and ! equal, then increment
            
            //  Rprintf( "same=%d  len=%d  count=%d\n", same, len, count );                
            
            //  start a new run
            ifirst      = i;
            valfirst    = pt[i] ;
            }
        }
    
    R_Free( interior );
    
    SEXP    out = PROTECT( Rf_allocVector(INTSXP,1) );
    *INTEGER(out)   =  count ;
    UNPROTECT(1) ;
    
    return( out ) ;
    }
    