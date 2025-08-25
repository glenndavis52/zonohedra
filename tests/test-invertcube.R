
library( zonohedra )
library( arrangements )

options( width=160 )   


#   gens    # of generators, based on polarzonohedron()
#   mults   # of the generators to turn into multiple groups, selected at random
#   reps    # to put in each group
#   loops   of loops to add, selected at random

randomzono <- function( gens, mults=as.integer(gens/2), reps=10, loops=4 )
    {
    set.seed(0)

    mat     = getmatrix( polarzonohedron(gens) ) #;  print( mat )

    #  choose mults points at random, and replicate them with multiples
    idxmult = arrangements::combinations( ncol(mat), mults, nsample=1 )

    idx     = unlist( lapply( 1:ncol(mat), function(i) { if( i %in% idxmult ) rep(i,reps) else i } ) )

    mat     = mat[ , idx ] #;  print( mat )

    #   now apply random scaling
    mat = matrix( runif(ncol(mat)), nrow=nrow(mat), ncol=ncol(mat), byrow=TRUE ) *  mat   #; print(mat)

    #   add loops
    idx = integer( ncol(mat) )
    #print( combinations(length(idx),loops) )

    idx[ arrangements::combinations(length(idx),loops,nsample=1) ] = 1L
    idx = (1:length(idx)) + cumsum(idx)

    #   return(idx)

    out     = matrix( 0, nrow(mat), ncol=idx[length(idx)] )

    out[ , idx ]  = mat #; print(out)

    out = zonohedron( out )

    return( out )
    }


#   returns  m x count matrix with random 2-trans cube points in the rows

random2trans    <- function( m, count )
    {
    set.seed(0)

    out = matrix( 0, count, m )

    trans = integer( count )

    for( i in 1:count )
        {
        pcube = numeric(m)

        #   choose a random pair - for point in pgram
        pair    = arrangements::combinations(m,2,nsample=1)

        pcube[ pair ]   = runif(2)

        if( pair[1]+1 <= pair[2]-1 )
            pcube[ (pair[1]+1) : (pair[2]-1) ]  = 1

        r   = runif(1)
        if( r <= 0.1 )
            pcube[ pair[1] ] = 0    # push to edge
        else if( r <= 0.2 )
            pcube[ pair[2] ] = 0    # push to other edge
        else if( r <= 0.3 )
            {
            #   push to vertex, this might make a vertex with only 0 transitions
            pcube[ pair[1] ] = 0  
            pcube[ pair[2] ] = 0              
            }

        if( runif(1) <= 0.5 )
            #   complement
            pcube   = 1 - pcube

        out[i, ]    = pcube

        trans[i]    = zonohedra:::transitioncount( pcube )

        if( 2L < trans[i] )
            {
            mess = sprintf( "random2trans().  trans=%d > 2.\n", trans[i] )
            cat( mess, file=stderr() )
            return(NULL)
            }
        }

    # out = cbind( out, trans )

    return( out )
    }

    
#   return TRUE or FALSE    
testinvert  <- function( tol=5.e-12 )
    {
    gens    = 17L
    
    zono = randomzono( gens )
    if( is.null(zono) ) return(FALSE)
    
    count = 1000

    pcube = random2trans( gens, count )
    if( is.null(pcube) ) return(FALSE)
    
    pcubelift   = zonohedra:::invertcubepoints( zono, pcube )
    
    # verify # of transitions
    for( i in 1:count )
        {
        trans = zonohedra:::transitioncount( pcubelift[i, ] )
        if( 2L < trans )
            {
            mess = sprintf( "testinvert(). i=%d  trans=%d > 2.\n", i, trans )
            cat( mess, file=stderr() )
            print( pcube[i, ] )
            print( pcubelift[i, ] )
            return(FALSE)
            }
        }
        
    # project and compare to original
    pcube_back = zonohedra:::projectcubepoints( zono, pcubelift )
    
    if( tol < max( abs(pcube_back - pcube) ) )
        {
        mess = sprintf( "testinvert(). max( abs(pcube_back - pcube) ) = %g > %g.\n", 
                        max( abs(pcube_back - pcube) ), tol  )
        cat( mess, file=stderr() )
        return( FALSE )
        }
    
    return(TRUE)
    }
    
    
if( ! testinvert() )  stop( "testinvert() failed !" )
    
cat( "Passed all cube inversion tests !\n", file=stderr() )       