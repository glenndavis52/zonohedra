require( zonohedra, quiet=TRUE )

options( width=160 )


g.zonolist  = list()

#g.zonolist[[1]] = zonohedron( classics.genlist[[3]] )
#names( g.zonolist )[1]  = attr( g.zonolist[[1]], "fullname" )

#g.zonolist[[2]] = zonohedron( classics.genlist[[6]] )
#names( g.zonolist )[2]  = attr( g.zonolist[[2]], "fullname" )

#g.zonolist[[3]] = zonohedron( classics.genlist[[13]] )
#names( g.zonolist )[3]  = attr( g.zonolist[[3]], "fullname" )

#g.zonolist[[4]] = zonohedron( colorimetry.genlist[[1]] )
#names( g.zonolist )[3]  =  attr( g.zonolist[[4]], "fullname" )

#g.zonolist[[1]] = zonohedron( colorimetry.genlist[[2]] )
#names( g.zonolist )[1]  =  attr( g.zonolist[[1]], "fullname" )

#g.zonolist[[4]] = spherize( g.zonolist[[1]] )
#names( g.zonolist )[4]  = paste( names[1], "after spherized", sep=' ', collapse='' )



#print( str(g.zonolist) )







#   returns a data frame with these columns
#       idxpair
#       alpha
#       direction

randomDirs <- function( matgen, base, level=0 )
    {
    out = randomAtLevel( matgen, level )

    out$direction   = out$XYZ  - matrix( base, nrow=nrow(out$XYZ), ncol=ncol(out$XYZ), byrow=TRUE )
    rownames( out$direction )   = 1:nrow(out$direction)

    out$XYZ     = NULL

    #rownames(out)   = 1:nrow(out)

    return( out )
    }



randomAtLevel <- function( matgen, level=0 )
    {
    set.seed(0)

    n   = ncol(matgen)

    level   = as.integer(level)

    ok  = (0 <= level) && (level <= n-2)
    if( ! ok )  return(NULL)

    #white   = rowSums( matgen )

    #   out = matrix( 0, 2*n, 3 )

    idxpair = matrix( NA_integer_, n, 2 )
    alpha   = matrix( NA_real_, n, 2 )
    XYZ     = matrix( NA_real_, n, 3 )

    for( k in 1:n )
        {
        shift   = k-1L

        #if( shift == 0 )
        #    perm    = 1:n
        #else
        #    perm    = c( (shift+1):n,1:shift )   # c( (n-shift+1):n, 1:(n-shift) )

        #print( perm )

        j1 = 1L
        j2 = j1 + level + 1L

        idxpair[k, ]    = ( ( c(j1,j2) + shift-1L ) %% n ) + 1L     #     perm[ c(j1,j2) ]
        alpha[k, ]      = runif(2)

        #pcube   = numeric(n)
        #pcube[ c(j1,j2) ]   = alpha[k, ]
        #if( 0 < level ) pcube[ (j1+1):(j2-1) ]  = 1
        #pcube[perm]   = pcube   # rotate cyclic

        pcube   = zonohedra:::pcubefromdata( idxpair[k, ], alpha[k, ], n )
        #print( pcube )

        # and multiply
        XYZ[k, ] = as.double( matgen %*% pcube )

        #   now add the antipodal data
        #idxpair[k+n, ]  = perm[ c(j2,j1) ]
        #alpha[k+n, ]    = 1 - alpha[k, ]
        #XYZ[k+n, ]      = white - XYZ[k, ]
        }

    out = data.frame( row.names=1:(n) )
    out$idxpair = idxpair
    out$alpha   = alpha
    out$XYZ     = XYZ

    return( out )
    }



randomCubeInterior <- function( n, reset=TRUE )
    {
    if( reset ) set.seed(0)

    eps = 5.e-4

    out = matrix( runif(n*3,min=eps,max=1-eps), n, 3 )

    return( out )
    }


randomCubeBoundary <- function( n )
    {
    set.seed(0)

    coord   = as.integer( runif(n,min=1,max=4) )

    value   = ( sign( runif(n,min=-1,max=1) ) + 1 ) / 2

    out = randomCubeInterior( n, FALSE )

    out[ cbind(1:n,coord) ]    = value

    return( out )
    }






testRaytrace <- function( tol=5.e-6 )
    {
    set.seed(0)

    g.zonolist[[4]] = zonohedron( colorimetry.genlist[['xyz1931.1nm']] )
    names( g.zonolist )[4]  =  attr( g.zonolist[[4]], "fullname" )


    for( k in 1:length(g.zonolist) )
        {
        zono    = g.zonolist[[k]]

        if( is.null(zono) ) next

        cat( "\n---------  testRaytrace() ", k, ")  ", names(g.zonolist)[k], '  --------------\n', sep='' )
        flush.console()

        basepoints  = 3
        #rays        = 10

        whitept = 2 * getcenter( zono )

        matsimp = getsimplified( getmatroid(zono) )

        matgen  = getmatrix( matsimp )

        m   = ncol(matgen)

        idxfromgnd  = zonohedra:::idxfromgroundfun( getground( matsimp ) )

        levelvec = c( 0, 1, as.integer(m/2), m-4 )

        levelvec = c( 1, as.integer(m/4), as.integer(m/2), as.integer(2*m/3), m-4 )
        #   levelvec = c( 1 )

        timetrace   = numeric(0)

        for( i in 1:basepoints )
            {
            s   = i / (basepoints+1)

            base = s * whitept

            for( level in levelvec )
                {
                dftest  = randomDirs( matgen, base, level )

                df  = raytrace2trans( zono, base, dftest$direction )     #; print( str(df) )

                #   check for failures
                if( any( is.na(df$tmax) ) )
                    {
                    cat( "testRaytrace(). base=", base, "  level=", level, '\n' )
                    cat( "   tmax failures=", sum( is.na(df$tmax) ), '\n' )
                    return(FALSE)
                    }

                #   compare idxpair values
                idxpair     = idxfromgnd[ df$gndpair ]

                if( any( dftest$idxpair != idxpair ) )
                    {
                    cat( "testRaytrace(). base=", base, "  level=", level, '\n' )
                    cat( "   idxpair failures=", sum( dftest$idxpair != idxpair ), '\n' )
                    return(FALSE)
                    }

                #   compare alpha values
                delta   = dftest$alpha - df$alpha

                good    = abs( delta ) <= tol
                if( any( ! good ) )
                    {
                    cat( "testRaytrace(). base=", base, "  level=", level, '\n' )
                    cat( "    alpha failures=", sum( ! good ), "  tol=", tol, '\n' )
                    cat( "    range(delta) = ", range(delta), '\n' )
                    return(FALSE)
                    }

                delta   = df$tmax - 1
                if( any( tol < abs(delta) ) )
                    {
                    count   = sum( tol < abs(delta) )
                    mess    = sprintf( "%d of the %d tmax values are outside tol=%g.\n",
                                            count, length(df$tmax), tol )
                    cat( mess )
                    return(FALSE)
                    }

                timetrace   = c( timetrace, df$timetrace )

                if( FALSE )
                    {
                    cat( '\n' )
                    mess = sprintf( "    fivenum summary of iters  --  total rays=%d  --\n", length(df$iters) )
                    cat( mess )
                    print( fivenum(df$iters) )
                    cat( "    mean(iters) =", mean(df$iters), '\n' )
                    flush.console()
                    }
                }
            }

        cat( '\n' )
        mess = sprintf( "fivenum summary of timetrace  --  total rays=%d  --\n", length(timetrace) )
        cat( mess )
        print( fivenum(timetrace) )

        cat( "    mean(timetrace) =", mean(timetrace), '\n' )
        }

    return( TRUE )
    }


testCubeRays <- function( basepoints=50, rays=200,  tol=5.e-13 )
    {
    cat( "\n----   testCubeRays()  -------\n" ) ; flush.console()

    cube    = zonohedron( classics.genlist[["C"]] )

    #   make basepoints
    basemat = randomCubeInterior( basepoints )

    #   make boundary points
    boundarymat     = randomCubeBoundary( rays )

    df  = NULL
    for( k in 1:basepoints )
        {
        base    = basemat[k, ]

        dirmat  = boundarymat - matrix( base, nrow(boundarymat), length(base), byrow=TRUE )

        dfsub   = raytrace2trans( cube, base, dirmat )

        #   check for any errors
        if( any( is.na(dfsub$tmax) ) )
            {
            count   = sum( is.na(dfsub$tmax) )
            mess    = sprintf( "For basepoint %g,%g,%g,  there were %d failures, out of %d.\n",
                                    base[1], base[2], base[3], count, length(dfsub$tmax) )
            cat( mess )
            return(FALSE)
            }


        if( tol < max( abs( boundarymat - dfsub$point ) ) )
            {
            count   = sum( tol < abs( boundarymat - dfsub$point ) )
            mess    = sprintf( "For basepoint %g,%g,%g,  %d of the %d boundary points are outside tol=%g.\n",
                                    base[1], base[2], base[3], count, nrow(boundarymat), tol )
            cat( mess )
            return(FALSE)
            }

        df  = rbind( df, dfsub )
        }

    delta   = df$tmax - 1
    if( any( tol < abs(delta) ) )
        {
        count   = sum( tol < abs(delta) )
        mess    = sprintf( "%d of the %d tmax values are outside tol=%g.\n",
                                count, length(df$tmax), tol )
        cat( mess )
        return(FALSE)
        }

    #print( df )


    mess = sprintf( "    fivenum summary of timetrace  --  total rays=%d  --\n", length(df$timetrace) )
    cat( mess )
    print( fivenum(df$timetrace) )
    cat( "    mean(timetrace) =", mean(df$timetrace), '\n' )

    cat( '\n' )
    mess = sprintf( "    fivenum summary of iters  --  total rays=%d  --\n", length(df$iters) )
    cat( mess )
    print( fivenum(df$iters) )
    cat( "    mean(iters) =", mean(df$iters), '\n' )

    return(TRUE)
    }


testPoles <- function( basepoints=50, tol=5.e-9 )
    {
    cat( "\n----   testPoles()  -------\n" ) ; flush.console()

    set.seed(0)

    zono    = zonohedron( classics.genlist[["BD"]] )

    matgen  = getmatrix(zono)

    n   = ncol( matgen )

    #   make random points in n-cube
    alpha   = matrix(  runif( n*basepoints,min=0.001, max=0.999 ), nrow=n, ncol=basepoints )

    #   make basepoints
    basemat = t( matgen %*% alpha )     #; print( basemat )

    #   make rays through "white"
    white   = 2 * getcenter(zono)
    #direction   = matrix( white, nrow=nrow(basemat), ncol=3, byrow=TRUE ) - basemat

    #   make rays through 0
    #direction   = rbind( direction, -basemat )

    for( k in 1:nrow(basemat) )
        {
        base    = basemat[k, ]

        direction   = rbind( white - base, -base )

        df  = raytrace2trans( zono, base, direction, tol=tol )  #; print( df )

        #   gndpair must be all NA
        if( ! all( is.na(df$gndpair) ) )
            {
            cat( "testPoles().  gndpair is invalid.\n" )
            print( df$gndpair )
            return(FALSE)
            }

        delta   = df$tmax - 1   #; print( range(delta) )

        if( any( tol < abs(delta) ) )
            {
            count   = sum( tol < abs(delta) )
            mess    = sprintf( "testPoles().  %d of the %d tmax values are outside tol=%g.\n",
                                    count, length(df$tmax), tol )
            cat( mess )
            return(FALSE)
            }
        }

    return( TRUE )
    }

#   1600 munsell chips

testMunsell <- function( path=NULL, tol=5.e-6  )
    {
    cat( "\n----   testMunsell()  -------\n" ) ; flush.console()

    if( is.null(path) ) path = system.file( "extdata/Munsell1600.txt", package='zonohedra' )

    timestart   = zonohedra:::gettime()

    zono = zonohedron( colorimetry.genlist[['xyz1931.1nm']] )

    #   calibrate
    sums    = rowSums( getmatrix(zono) )  #; print(sums)

    #   use a diagonal scaling to match Illuminant E
    zono = lintransform( zono, diag(100/sums) )

    XYZ = as.matrix( read.table( path ) ) # ; print( str(XYZ) )

    #XYZ = XYZ[ 1492:1541, ]
    #print( XYZ / XYZ[ ,2] )

    base    = getcenter(zono) #; print(2*base)

    direction   = XYZ - matrix(base,nrow=nrow(XYZ),ncol=ncol(XYZ),byrow=TRUE)

    df  = raytrace2trans( zono, base, direction, tol=tol )  #; print( df )

    #   check for any errors
    if( any( is.na(df$tmax) ) )
        {
        count   = sum( is.na(df$tmax) )
        mess    = sprintf( "For basepoint %g,%g,%g,  there were %d failures, out of %d.\n",
                                base[1], base[2], base[3], count, length(df$tmax) )
        cat( mess )
        return(FALSE)
        }

    cat( "time elapsed = ", zonohedra:::gettime() - timestart, "sec  \n" )
    cat( '\n' )

    mess = sprintf( "    fivenum summary of timetrace  --  total rays=%d  --\n", length(df$timetrace) )
    cat( mess )
    print( fivenum(df$timetrace) )
    cat( "    mean(timetrace) =", mean(df$timetrace), "   sum(timetrace) =", sum(df$timetrace),'\n' )

    cat( '\n' )
    mess = sprintf( "    fivenum summary of iters  --  total rays=%d  --\n", length(df$iters) )
    cat( mess )
    print( fivenum(df$iters) )
    cat( "    mean(iters) =", mean(df$iters), '\n' )


    return(TRUE)
    }


testSections <- function( genlimit=100, tol=5.e-9 )
    {
    set.seed(0)

    g.zonolist[[1]] = zonohedron( classics.genlist[[3]] )

    g.zonolist[[2]] = zonohedron( classics.genlist[[6]] )

    g.zonolist[[3]] = zonohedron( classics.genlist[[8]] )

    g.zonolist[[4]] = zonohedron( colorimetry.genlist[['xyz1931.1nm']] )

    for( k in 1:length(g.zonolist) )
        {
        zono    = g.zonolist[[k]]

        if( is.null(zono) ) next

        cat( '\n\n' )
        cat( "---------  testSections() ", k, ")   ", attr(zono,"fullname"), ' --------------\n', sep='' )
        flush.console()

        directions  = 20

        if( genlimit <= ncol( getmatrix(zono) ) )   directions = 2

        sections    = 5

        #   make uniformly random directions
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

            res = section2trans( zono, normal, beta, invert=TRUE )

            count   = length(res)

            #   count the number of points
            for( j in seq_len(count) )
                points  = points + nrow(res[[j]]$point)

            #  if( genlimit <= ncol( getmatrix(zono) ) )    count=0  # too big for verification and inversion

            #   do some verification and inversion
            for( j in seq_len(count) )
                {
                df      = res[[j]]

                point   = df$point

                m   = nrow(point)

                #   cat( "    direction=", i, "  beta=", beta[j], "  points=", m, '\n' )

                #   all points must be not NA
                mask        = is.finite( point )
                ok  = all( mask )
                if( ! ok )
                    {
                    cat( "testSections(). k=", k, "  i=", i, "  j=", j, "  normal=", normal, "\n" )
                    cat( "   section finite failures=", sum( ! mask ), '\n' )
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

                #   check inversion accuracy
                delta   =  pcube  %*%  t(zono$matroid$matrix) -  point   #; print(delta)
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


if( ! testPoles() )     stop( "testPoles() failed !" )

if( ! testCubeRays() )  stop( "testCubeRays() failed !" )

if( ! testRaytrace() )  stop( "testRaytrace() failed !" )

if( ! testMunsell() )  stop( "testMunsell() failed !" )

if( ! testSections() )  stop( "testSections() failed !" )


cat( "\nPassed all 2-transition tests !\n", file=stderr() )
