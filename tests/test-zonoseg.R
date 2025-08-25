require( zonohedra, quiet=TRUE )
require( arrangements, quiet=TRUE )
options( width=144 )   

testZonosegInterior <- function( nonneg, tol=1.e-12 )
    {
    set.seed(0)
    
    for( k in 1:20 )
        {
        #   make a random zonoseg 
        points  = 100
        loops   = 10
        
        gen =   rnorm(points)
        
        if( nonneg )    gen = abs(gen)
        
        idxloop = arrangements::combinations( points, loops, nsample=1 )
        gen[ idxloop ] = 0
        
        zono    = zonoseg( gen )
        if( is.null(zono) ) return(FALSE)
        
        #   make random points in the segment
        n   = 100
        seg = getsegment(zono)
        z   = runif( n, min=seg[1], max=seg[2] )
        
        pcube  = invert( zono, z )$pcube
        
        zback   = pcube %*% gen
        
        delta   = zback - z     #;  cat( "range(delta) = ", range(delta), '\n' )
        
        ok  = all( abs(delta) <= tol )
        if( ! ok )
            {
            cat( "testZonosegInterior().  ERR.  range(delta) = ", range(delta), '\n' )
            return(FALSE)
            }
        }

    return(TRUE)
    }
    
    
testZonosegExterior <- function( )
    {
    set.seed(0)
    
    for( k in 1:20 )
        {
        #   make a random zonoseg 
        points  = 100
        loops   = 10
        
        gen =   rnorm(points)
        
        #   if( nonneg )    gen = abs(gen)
        
        idxloop = arrangements::combinations( points, loops, nsample=1 )
        gen[ idxloop ] = 0
        
        zono    = zonoseg( gen )
        if( is.null(zono) ) return(FALSE)
        
        #   make random points on either side of the segment
        n   = 50
        seg = getsegment(zono)
        z   = c( runif( n, min=seg[1]-1, max=seg[1] ) , runif( n, min=seg[2], max=seg[2]+1 ) )
        
        pcube  = invert( zono, z )$pcube

        ok  = all( is.na(pcube) )
        if( ! ok )
            {
            cat( "testZonosegExterior(). ERR.  some cube values are finite.\n" )
            return(FALSE)
            }
        }

    return(TRUE)
    }
    
        
testZonosegBoundary <- function( tol=1.e-12 )
    {
    set.seed(0)
    
    for( k in 1:20 )
        {
        #   make a random zonoseg 
        points  = 100
        loops   = 10
        
        gen =   rnorm(points)

        idxloop = arrangements::combinations( points, loops, nsample=1 )
        gen[ idxloop ] = 0
        
        zono    = zonoseg( gen )
        if( is.null(zono) ) return(FALSE)
        
        #   test both endpoints
        z   = getsegment(zono)

        pcube  = invert( zono, z )$pcube
        
        zback   = pcube %*% gen
        
        delta   = zback - z     #;  cat( "range(delta) = ", range(delta), '\n' )
        
        ok  = all( abs(delta) <= tol )
        if( ! ok )
            {
            cat( "testZonosegBoundary().  ERR.  range(delta) = ", range(delta), '\n' )
            return(FALSE)
            }
        }

    return(TRUE)
    }
    
if( ! testZonosegInterior(TRUE) )  stop( "testZonosegInterior(TRUE) failed !" )
    
if( ! testZonosegInterior(FALSE) )  stop( "testZonosegInterior(FALSE) failed !" )

if( ! testZonosegExterior() )  stop( "testZonosegExterior() failed !" )

if( ! testZonosegBoundary(tol=0) )  stop( "testZonosegBoundary() failed !" )



cat( "Passed all zonoseg tests !\n", file=stderr() )    
