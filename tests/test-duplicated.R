
library( zonohedra )
library( microbenchmark )

options( width=120 )

testDuplicated <- function()
    {
    set.seed(0)
    
    start   = 1000
    reps    = 10
    A       = matrix( rnorm(3*start), 3, start*reps )

    cat( "A = \n" )
    print( str(A) )


    cat( "---------- duplicated  (std::set  or  std::unordered_set) ---------------\n" )
    res     = base::duplicated.matrix(A,MARGIN=2)
    print( str(res) )
    print( res[ (ncol(A)-10):ncol(A) ] )

    res     = zonohedra::duplicated.matrix(A,MARGIN=2)
    print( str(res) )
    print( res[ (ncol(A)-10):ncol(A) ] )

    print( microbenchmark( base::duplicated.matrix(A,MARGIN=2)  ) )
    print( microbenchmark( zonohedra::duplicated.matrix(A,MARGIN=2) ) )

    return(TRUE)
    }
    
    
testGrpDuplicated <- function( start=1000, reps=10 )
    {
    cat( "\n" )
    
    A       = matrix( rnorm(3*start), 3, start*reps )
    
    cat( "---------- group duplicated (std::map  or  std::unordered_map)  ---------------\n" )
    grouptarget = rep( 1:start, reps )
    group       = zonohedra::grpDuplicated(A,MARGIN=2)
    if( ! all(group==grouptarget) )
        {
        cat( "group incorrect:\n" )
        print(group)
        return(FALSE)
        }
        
    #print( str(res) )
    #print( res[ (ncol(A)-10):ncol(A) ] )
    print( microbenchmark( zonohedra::grpDuplicated(A,MARGIN=2) )  )

    return(TRUE)
    }
    
    

# if( ! testDuplicated() )  stop( "testDuplicated() FAILED !  ERROR" )
 
if( ! testGrpDuplicated() )  stop( "testGrpDuplicated() FAILED !  ERROR" )
 
 
cat( "Passed all duplication tests !\n", file=stderr() )        
