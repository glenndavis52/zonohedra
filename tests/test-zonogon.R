require( zonohedra, quiet=TRUE )

options( width=144 )

set.seed(0)


g.zonolist  = list()

g.zonolist[[1]] = zonogon( matrix( 1:4, nrow=2, ncol=2 ) )
names( g.zonolist )[1]  = "parallelogram"

g.zonolist[[2]] = symmetrize( g.zonolist[[1]] )
names( g.zonolist )[2]  = "symmparallelogram"

g.zonolist[[3]] = polarzonogon( 17 )
names( g.zonolist )[3]  = "polar"

g.zonolist[[4]] = polarzonogon( 17, 9 )
names( g.zonolist )[4]  = "halfpolar"

mat = matrix( rnorm(2*50), nrow=2 )
g.zonolist[[5]] = zonogon( mat )
names( g.zonolist )[5]  = "random"

mat = rbind( rnorm(25), runif(25) )
g.zonolist[[6]] = zonogon( mat )
names( g.zonolist )[6]  = "halfrandom"


#print( str(g.zonolist) )


testInsideAndInversion <- function( tol=5.e-12 )  # before May 2022 it was smaller 1.e-14, but this failed for the "random" case
    {
    set.seed(0)

    for( k in 1:length(g.zonolist) )
        {
        cat( '\n' )
        cat( "---------  testInsideAndInversion() ", k, ")  ", names(g.zonolist)[k], '--------------\n' )

        zono    = g.zonolist[[k]]

        basepoints  = 9
        rays        = 10

        whitept = 2 * zono$center

        for( i in 1:(basepoints-1) )
            {
            s   = i / (basepoints)

            base = s * whitept

            direction   = matrix( rnorm(2*rays), nrow=rays )

            df  = raytrace( zono, base, direction )     #; print( str(df) )
            point    = df$point
            
            #print( df )
            
            #   check that distance of these boundary points is very small
            distance    = inside( zono, point )$distance #;    print(distance)

            mask    = abs(distance) <= tol

            ok  = all( mask, na.rm=TRUE )
            if( ! ok )
                {
                cat( "testInsideAndInversion(). k=", k, "  i=", i, "  tol=", tol, "\n" )
                cat( "   raytrace boundary failures=", sum( ! mask, na.rm=TRUE ),  "  max(abs(distance))=", max(abs(distance)), '\n' )
                return(FALSE)
                }
                
            #   test inversion of these boundary points
            #print( point )
            pcube   = invert( zono, point, tol=tol )$pcube
            if( is.null(pcube) )    
                {
                cat( "testInsideAndInversion(). k=", k, "  i=", i, "  tol=", tol, "\n" )
                cat( "pcube is NULL, for boundary points.\n" )                
                return(FALSE)
                }
                
            #print( pcube )
                
            delta   =  pcube  %*%  t(zono$matroid$matrix) -  point   #; print(delta)
            ok  = all( abs(delta) < tol )   #; print(ok)
            if( ! ok )
                {
                cat( "testInsideAndInversion(). k=", k, "  i=", i, "\n" )
                cat( "   inversion failures for boundary points=", sum( ! (abs(delta) < tol), na.rm=TRUE ), "   max=", max(abs(delta)), '\n' )
                print( pcube )
                return(FALSE)
                }
                

            #   on each segment from base to boundary,
            #   pick a point at random and make sure the new points are all inside
            lambda  = runif( nrow(point) )
            pmat    = lambda * matrix( base, nrow=nrow(point), ncol=2, byrow=TRUE )  +  (1-lambda) * point
            inside  = inside( zono, pmat )$inside #; print(inside)

            ok  = all( inside, na.rm=TRUE )
            if( ! ok )
                {
                cat( "testInsideAndInversion(). k=", k, "  i=", i, "\n" )
                cat( "   inside failures=", sum( ! inside, na.rm=TRUE ), '\n' )
                return(FALSE)
                }
                
            #   test inversion of these inside points
            pcube   = invert( zono, pmat, tol=tol )$pcube
            if( is.null(pcube) )    
                {
                cat( "testInsideAndInversion(). k=", k, "  i=", i, "  tol=", tol, "\n" )
                cat( "pcube is NULL, for inside points.\n" )                
                return(FALSE)
                }
                
            delta   =  pcube  %*%  t(zono$matroid$matrix) -  pmat   #; print(delta)
            ok  = all( abs(delta) < tol )   #; print(ok)
            if( ! ok )
                {
                cat( "testInsideAndInversion(). k=", k, "  i=", i, "\n" )
                cat( "   inversion failures on inside =", sum( ! (abs(delta) < tol), na.rm=TRUE ), "   max=", max(abs(delta)), '\n' )
                print( pcube )
                return(FALSE)
                }

            #   move a little past the boundary and make sure the new points are all outside
            test    = -0.1 * matrix( base, nrow=nrow(point), ncol=2, byrow=TRUE )  +  1.1 * point
            outside = ! inside( zono, test )$inside  #; print(outside)

            ok  = all( outside, na.rm=TRUE )
            if( ! ok )
                {
                cat( "testInsideAndInversion(). k=", k, "  i=", i, "\n" )
                cat( "   outside failures=", sum( ! outside, na.rm=TRUE ), '\n' )
                return(FALSE)
                }
            }
        }

    return(TRUE)
    }


testSections <- function( tol = 1.e-14)
    {
    set.seed(0)

    for( k in 1:length(g.zonolist) )
        {
        cat( '\n' )
        cat( "---------  testSections() ", k, ")  ", names(g.zonolist)[k], '--------------\n' )

        zono    = g.zonolist[[k]]

        directions  = 100

        direction   = matrix( rnorm(2*directions), nrow=directions )

        for( i in 1:directions )
            {
            normal  = direction[i, ]

            #   find max and min of the functional over the zonogon
            df  = support( zono, rbind(normal,-normal) )
            valmax  =  df$value[1]
            valmin  = -df$value[2]

            #   cat( "normal=", normal,  "  valmin=", valmin, "valmax=", valmax, '\n' )

            #   generate 9 sections between min and max
            s   = (1:9)/10
            beta    = (1-s)*valmin + s*valmax
            df  = section( zono, normal, beta )
            boundary1   = df$boundary1
            boundary2   = df$boundary2
            
            #   all points must be not NA
            point       = rbind(boundary1,boundary2)
            mask        = is.finite( point )
            ok  = all( mask )
            if( ! ok )
                {
                cat( "testSections(). k=", k, "  i=", i, "  normal=", normal, "\n" )
                cat( "   section failures=", sum( ! mask ), '\n' )
                return(FALSE)
                }

            #   check that distance of these boundary points is very small
            distance    = inside( zono, point )$distance     #;    print(distance)

            mask    = abs(distance) <= tol

            ok  = all( mask )
            if( ! ok )
                {
                cat( "testSections(). k=", k, "  i=", i, "  normal=", normal, "  tol=", tol, "\n" )
                cat( "   boundary failures=", sum( ! mask ), '\n' )
                return(FALSE)
                }

            #   make a few interpolation points, they must all be inside
            for( s in c(0.01,0.5,0.99) )
                {
                test    = (1-s)*boundary1 + s*boundary2
                inside  = inside( zono, test )$inside
                
                ok  = all( inside  )
                if( ! ok )
                    {
                    cat( "testSections(). k=", k, "  i=", i, "  normal=", normal, "\n" )
                    cat( "   inside failures=", sum( ! inside ), '\n' )
                    return(FALSE)
                    }
                }
                
            #   make a few extrapolation points, they must all be outside
            for( s in c(-1,-0.01,1.01,2) )
                {
                test    = (1-s)*boundary1 + s*boundary2
                outside = ! inside( zono, test )$inside
                
                ok  = all( outside  )
                if( ! ok )
                    {
                    cat( "testSections(). k=", k, "  i=", i, "  normal=", normal, "\n" )
                    cat( "   outside failures=", sum( ! outside ), '\n' )
                    return(FALSE)
                    }
                }
            }
        }

    return(TRUE)
    }

testSums <- function()
    {
    thesum = g.zonolist[[1]] %+%  g.zonolist[[2]] %+%  g.zonolist[[3]] %+%  g.zonolist[[4]] %+%  g.zonolist[[5]] %+%  g.zonolist[[6]]

    if( is.null(thesum) )   return(FALSE)
    
    print(thesum)
    
    return(TRUE)
    }
    
if( ! testSums() )  stop( "testSums() failed !" )

if( ! testInsideAndInversion() )  stop( "testInsideAndInversion() failed !" )
    
if( ! testSections() )  stop( "testSections() failed !" )

cat( "\nPassed all zonogon tests !\n", file=stderr() )

warnings()