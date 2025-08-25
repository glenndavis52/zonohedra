
#include <stdbool.h>

//  vec     vector of doubles
//  m       length of vec

//  return the mean of the values in vec[]
//  Exception: If vec contains a single integer (possibly with repeats),
//  then return that integer.

static double
collapsedValue( const double *vec, int m )
    {
    bool    testint = true ;
    bool    haveint = false ;   //  a single integer has been found
    int     integer = 0;            //  the integer that was found in vec[]
        
    double  sum = 0 ;
    
    for( int k=0 ; k<m ; k++ )
        {
        sum += vec[k];
        
        if( testint && ((int)vec[k])==vec[k] )
            {
            //  vec[k] is an integer
            if( haveint && vec[k]!=integer )
                {
                //  found a 2nd integer, so cancel everything
                haveint = false ;
                testint = false ;
                }
            else
                {
                //  found the 1st integer, or a repeat
                haveint = true ;
                integer = vec[k];
                }
            }
        }
    
    if( haveint )
        return( integer );
    else
        return( sum / m );
    }
    

//  Exception: if the group contains a single integer (possibly with repeats),
//  then each value in the group is replaced by that integer.


//  vec     vector of doubles sorted in ascending order
//  n       length of vec
//  eps     positive small constant

//  A _group_ is a maximal set of consecutive elements in vec
//  with adjacent differences all <= eps

//  the function modifies vec[] by replacing each value in a group
//  by the mean of that group.
//  Exception: if the group contains a single integer (possibly with repeats),
//  then each value in the group is replaced by that integer.
//
//  returns true or false

bool
collapseGroups1D( double *vec, int n, double eps )
    {
    bool    ingroup = false ;
    int     ifirst = -1;    
    
    for( int i=1 ; i<n ; i++ )
        {
        double  delta = vec[i] - vec[i-1] ;
        
        if( delta <= eps )
            {
            if( ! ingroup )
                {
                //  start new group
                ingroup   = true ;                
                ifirst = i-1 ;
                }
            }
        else if( ingroup )
            {
            //  close group we just left
            double  val = collapsedValue( vec+ifirst, i-ifirst );
            ingroup = false ;                
            for( int k=ifirst ; k<i ; k++ )
                vec[k] = val ;
            }
        }
        
    if( ingroup )
        {
        //  close the last group
        double  val = collapsedValue( vec+ifirst, n-ifirst );        
        for( int k=ifirst ; k<n ; k++ )
            vec[k] = val ;
        }
        
    return true;
    }
    