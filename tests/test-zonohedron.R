require( zonohedra, quiet=TRUE )

options( width=160 )

set.seed(0)


g.zonolist  = list()

g.zonolist[[1]] = zonohedron( classics.genlist[[3]] )
names( g.zonolist )[1]  = attr( g.zonolist[[1]], "fullname" )

g.zonolist[[2]] = zonohedron( classics.genlist[[6]] )
names( g.zonolist )[2]  = attr( g.zonolist[[2]], "fullname" )

g.zonolist[[3]] = zonohedron( classics.genlist[[13]] )
names( g.zonolist )[3]  = attr( g.zonolist[[3]], "fullname" )

#g.zonolist[[4]] = zonohedron( colorimetry.genlist[[1]] )
#names( g.zonolist )[3]  =  attr( g.zonolist[[4]], "fullname" )

#   for the 1nm data, inside.zonohedron() takes too long !
g.zonolist[[4]] = zonohedron( colorimetry.genlist[[2]] )
names( g.zonolist )[4]  =  attr( g.zonolist[[4]], "fullname" )

#g.zonolist[[4]] = spherize( g.zonolist[[1]] )
#names( g.zonolist )[4]  = paste( names[1], "after spherized", sep=' ', collapse='' )



#print( str(g.zonolist) )


testRaytrace <- function( tol=5.e-13 )
    {
    set.seed(0)

    for( k in 1:length(g.zonolist) )     # )
        {
        cat( '\n' )
        cat( "---------  testRaytrace() ", k, ")  ", names(g.zonolist)[k], '--------------\n' )

        zono    = g.zonolist[[k]]

        basepoints  = 8
        rays        = 10

        whitept = 2 * zono$center

        for( i in 1:basepoints )
            {
            s   = i / (basepoints+1)

            base = s * whitept

            direction   = matrix( rnorm(3*rays), nrow=rays )

            df  = raytrace( zono, base, direction, invert=TRUE )     #; print( str(df) )

            point   = df$point

            #   check that distance of these boundary points is very small
            distance    = inside( zono, point )$distance #;    print(distance)

            mask    = abs(distance) <= tol

            ok  = all( mask, na.rm=TRUE )
            if( ! ok )
                {
                cat( "testRaytrace(). k=", k, "  i=", i, "  tol=", tol, "\n" )
                cat( "   boundary failures=", sum( ! mask, na.rm=TRUE ), '\n' )
                return(FALSE)
                }


            #   test inversion of these points on the boundary
            pcube   = df$pcube      # invertboundary( zono, point, tol=tol )$pcube
            if( is.null(pcube) )
                {
                cat( "testRaytrace(). k=", k, "  i=", i, "  tol=", tol, "\n" )
                cat( "    pcube is NULL, for boundary points.\n" )
                return(FALSE)
                }

            delta   =  pcube  %*%  t(zono$matroid$matrix) -  point   #; print(delta)
            ok  = all( abs(delta) < tol )   #; print(ok)
            if( ! ok )
                {
                cat( "testRaytrace(). k=", k, "  i=", i, "\n" )
                cat( "   inversion failures for boundary points=",
                        sum( ! (abs(delta) < tol), na.rm=TRUE ), "   max=", max(abs(delta)), '\n' )
                print( pcube )
                return(FALSE)
                }



            #   move a little from the boundary to the base and make sure the new points are all inside
            test    = 0.1 * matrix( base, nrow=nrow(point), ncol=3, byrow=TRUE )  +  0.9 * point
            inside  = inside( zono, test )$inside #; print(inside)

            ok  = all( inside, na.rm=TRUE )
            if( ! ok )
                {
                cat( "testRaytrace(). k=", k, "  i=", i, "\n" )
                cat( "   inside failures=", sum( ! inside, na.rm=TRUE ), '\n' )
                return(FALSE)
                }

            #   move a little past the boundary and make sure the new points are all outside
            test    = -0.1 * matrix( base, nrow=nrow(point), ncol=3, byrow=TRUE )  +  1.1 * point
            outside = ! inside( zono, test )$inside  #; print(outside)

            ok  = all( outside, na.rm=TRUE )
            if( ! ok )
                {
                cat( "testRaytrace(). k=", k, "  i=", i, "\n" )
                cat( "   outside failures=", sum( ! outside, na.rm=TRUE ), '\n' )
                return(FALSE)
                }
            }
        }

    return(TRUE)
    }


testSections <- function( genlimit=100, tol=5.e-9 )
    {
    set.seed(0)
    

    for( k in 1:length(g.zonolist) )
        {
        cat( '\n\n' )
        cat( "---------  testSections() ", k, ")   ", names(g.zonolist)[k], ' --------------\n', sep='' )
        flush.console()

        zono    = g.zonolist[[k]]

        directions  = 20
        
        if( genlimit <= ncol( getmatrix(zono) ) )   directions = 1

        sections    = 5

        direction   = matrix( rnorm(3*directions), nrow=directions )

        timeelapsed = 0
        timestart   = zonohedra:::gettime()

        points      = 0

        for( i in 1:directions )
            {
            normal  = direction[i, ]

            #   find max and min of the functional over the zonohedron
            df  = support( zono, rbind(normal,-normal) )
            valmax  =  df$value[1]
            valmin  = -df$value[2]

            #   cat( "normal=", normal,  "  valmin=", valmin, "valmax=", valmax, '\n' )

            #   generate sections between min and max
            s   = (1:sections)/(sections+1)

            beta    = (1-s)*valmin + s*valmax

            res = section( zono, normal, beta )

            count   = length(res)

            #   count the number of points
            for( j in seq_len(count) )
                points  = points + nrow(res[[j]]$point)

            #  if( genlimit <= ncol( getmatrix(zono) ) )    count=0  # too big for verification and inversion

            #   do some verification and inversion
            for( j in seq_len(count) )
                {
                section = res[[j]]$point

                m   = nrow(section)

                #   cat( "    direction=", i, "  beta=", beta[j], "  points=", m, '\n' )

                #   all points must be not NA
                mask        = is.finite( section )
                ok  = all( mask )
                if( ! ok )
                    {
                    cat( "testSections(). k=", k, "  i=", i, "  j=", j, "  normal=", normal, "\n" )
                    cat( "   section finite failures=", sum( ! mask ), '\n' )
                    return(FALSE)
                    }

                #   check that distance of these boundary points is very small
                distance    = inside( zono, section )$distance     #;    print(distance)

                mask    = abs(distance) <= tol * max( abs(valmax), abs(valmin) )

                ok  = all( mask )
                if( ! ok )
                    {
                    cat( "testSections(). k=", k, "  i=", i, "  j=", j, "  normal=", normal, "  tol=", tol, "\n" )
                    cat( "   boundary failures=", sum( ! mask ), '\n' )
                    cat( "   distance=", distance, '\n' )
                    return(FALSE)
                    }


                #   test inversion of these boundary points
                df  = invertboundary( zono, section, tol=tol )

                if( is.null(df) )
                    {
                    cat( "testSections(). k=", k, "  i=", i, "  tol=", tol, "\n" )
                    cat( "    pcube is NULL, for boundary points.\n" )
                    return(FALSE)
                    }

                pcube   = df$pcube
                
                if( any( is.na(pcube) ) )
                    {
                    cat( "testSections(). k=", k, "  i=", i, "\n" )
                    cat( "    normal =", normal, "  beta =", beta[j], '\n' )
                    mess    = sprintf( "    # of NAs in pcube is %d, out of %d values.\n", 
                                            sum(is.na(pcube)), length(pcube) )
                    cat( mess )
                    return(FALSE)
                    }                

                if( df$transitions[1] %% 2 == 1 )
                    {
                    cat( "    transitions[1]=", df$transitions[1], '\n' )
                    cat( "    pcube=", pcube[1, ], '\n' )
                    }

                delta   =  pcube  %*%  t(zono$matroid$matrix) -  section   #; print(delta)
                ok  = all( abs(delta) < tol )   #; print(ok)
                if( ! ok )
                    {
                    cat( "testSections(). k=", k, "  i=", i, "\n" )
                    cat( "   inversion failures for boundary points=",
                            sum( ! (abs(delta) < tol), na.rm=TRUE ), "   max=", max(abs(delta)), '\n' )
                    # print( pcube )
                    return(FALSE)
                    }

                }
            }

        timeelapsed = zonohedra:::gettime() - timestart

        sections    = directions * length(beta)

        cat( "    sections=", sections, "  time=", timeelapsed, "  ", timeelapsed/sections, " per section",
                  "  points=", points, "  ", timeelapsed/points, " per point.\n", sep='' )

        flush.console()
        }

    return(TRUE)
    }
    
    

library( arrangements )
library( zonohedra )


###     test raytracing of points in non-trivial facets of a zonohedron

#   thematroid      rank 3 matroid, does not have be simple
#   subgrnd         a subset of the ground set, must be contiguous (not verified)
#   pairs           the number of pairs to generate
#
#   returns a data.frame with pairs rows, and these columns:

buildmulti2trans <- function( thematroid, subgrnd, pairs )
    {
    set.seed(0)
    
    ground  = getground(thematroid)
    
    idxfromgnd  = zonohedra:::idxfromgroundfun( ground )

    idxpair = arrangements::combinations( idxfromgnd[subgrnd], 2, nsample=pairs ) #; print( idxpair )
    
    dim(idxpair)    = c(pairs,2)
    
    if( FALSE ) idxpair = rbind( idxpair, idxpair[ ,2:1,drop=F] )
    
    n   = length(ground) #; print( n )
    
    idxgnd      = ground[ idxpair ]
    dim(idxgnd) = dim(idxpair)

    pcube   = matrix( NA_real_, nrow=nrow(idxpair), ncol=n )
    colnames(pcube) = ground
    
    rank    = rep( NA_integer_, nrow(idxpair) )
    
    for( k in 1:nrow(idxpair) )
        {
        alpha   = runif(2)
        
        pcube[k, ]  = zonohedra:::pcubefromdata( idxpair[k, ], alpha, n )
        
        rank[k] = rank( thematroid, idxgnd[k,1] : idxgnd[k,2] )
        }
    
    out     = data.frame( row.names=1:nrow(idxpair) )
    out$idxpair = idxpair
    out$idxgnd  = idxgnd
    out$pcube   = pcube
    out$rank    = rank
    
    return(out)
    }
    
    
testRaytraceNT <- function( zono, pairs=200, tol1=1.e-4, tol2=1.e-2 )
    {
    cat( '\n' )
    cat( "---------  testRaytraceNT()  ",  attr(zono,"fullname"), '   --------------\n' )
    flush.console()

    
    ntcount = length(zono$zonogon)
    if( ntcount == 0 )  return(FALSE)
    
    #matsimple   = getsimplified(zono$matroid)
    
    thematroid  = getmatroid(zono)
    
    df  = buildmulti2trans( thematroid, thematroid$hyperplane[[1]], pairs )
    
    matgen  = getmatrix(thematroid)
    
    df$XYZ  = df$pcube  %*%  t(matgen)
    
    # print( str(df) )
    
    base        = getcenter(zono)
    direction   = df$XYZ - matrix(base,nrow=nrow(df$XYZ),ncol=3,byrow=T)
    
    dftest = raytrace( zono, base, direction, invert=TRUE ) 
    
    dftest$idxgnd   = df$idxgnd
    
    dftest$rank     = df$rank
    
    matgenorg   = getmatrix(zono)
    
    XYZ = dftest$pcube  %*%  t(matgenorg)
    
    dftest$deltaXYZ     = apply( abs(df$XYZ - XYZ), 1L, max )    # rowMeans( abs(df$XYZ - XYZ) )

    dftest$deltacube    = apply( abs(df$pcube - dftest$pcube), 1L, max )  # rowMeans( abs(df$pcube - dftest$pcube) )

    #   change deltacube entries with rank=1 to NA.  These are degenerate.
    mask    = df$rank<2
    dftest$deltacube[ mask ]  = NA_real_
    
    #print(dftest)
    
    mask    = dftest$transitions == 2
    ok  = all( mask )
    if( ! ok )
        {
        cat( "testraytraceNT(). FAIL.  transitions==2 failures=", sum( ! mask, na.rm=TRUE ),  
                    "   max(transitions) =", max( dftest$transitions ), '\n' )
        return(FALSE)
        }
        
    deltamax   = max( dftest$deltaXYZ, na.rm=T ) 
    if( tol1 < deltamax )
        {
        mess    = sprintf( "testraytraceNT().  FAIL.  max(deltaXYZ)=%g > %g.\n", deltamax, tol1 )
        cat( mess )
        return(FALSE)
        }
        
    deltamax    = max( dftest$deltacube, na.rm=T )
    if( tol2 < deltamax )
        {
        kmax    = which.max(dftest$deltacube)        
        mess    = sprintf( "testraytraceNT().  FAIL.  max(deltacube)=%g > %g.  argmax=%d,%d\n", 
                            deltamax, tol2, dftest$idxgnd[kmax,1], dftest$idxgnd[kmax,2] )
        cat( mess )
        return(FALSE)
        }
    
    return( TRUE )
    }    


if( ! testRaytrace() )          stop( "testRaytrace() failed !" )

#   for testing non-trivial facets, use the one facet in zonocol.1nm
if( ! testRaytraceNT( g.zonolist[[4]] ) )        stop( "testRaytraceNT() failed !" )

if( ! testSections( 400 ) )     stop( "testSections() failed !" )

cat( "\nPassed all zonohedron tests !\n", file=stderr() )
