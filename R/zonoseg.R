#
#   zonoseg is a 1-dimensional zonotope
#
#   implemented as a list with items:

#   matroid         the matroid for it, which includes the generating matrix, etc.

#   segment         the actual segment [ segment[1], segment[2] ]
#   center          midpoint of the segment
#   segment2trans   the subsegment coming from 2-transition points in the cube


#   zonoseg constructor
#
#   mat     a numeric matrix with 1 row
#   e0      threshold for a column vector to be considered 0


zonoseg <- function( mat, e0=0, ground=NULL )
    {
    if( ! is.matrix(mat) )
        {
        temp    = names(mat)
        mat = matrix( mat, nrow=1 )
        colnames(mat) = temp
        }

    ok  = is.matrix(mat)  &&  is.numeric(mat)  &&  nrow(mat)==1  &&  1<=ncol(mat)
    if( ! ok )
        {
        log_level( ERROR, "mat is invalid." )
        return(NULL)
        }

    matroid     = matroid( mat, e0=e0, ground=ground )
    if( is.null(matroid) )
        return(NULL)

    out = list()

    class( out )    = c( "zonoseg", "zonotope", class(out) )

    out$matroid = matroid

    #   shift is the sum of all negative entries
    #   and used later in invert()
    #out$shift   = sum( mat[ mat < 0 ] )

    generator   =  as.double(mat)

    df  = dfminmax( generator )
    out$segment = as.double(df$value)
    names(out$segment) = c("min","max")
    
    out$center  = mean( out$segment )

    df  = dftrans2( generator )
    out$segment2trans   = as.double(df$value)
    names(out$segment2trans) = c("min","max")

    return(out)
    }




print.zonoseg  <-  function( x, ... )
    {
    #   convert the 1-row matrix to a plain vector
    generator           = as.double( x$matroid$matrix )
    names(generator)    = as.character( x$matroid$ground )

    npoints = length(generator)

    idxfromground   = integer( max(x$matroid$ground) )
    idxfromground[ x$matroid$ground ] = 1:npoints

    idxloop = idxfromground[ x$matroid$loop ]

    maskloop    = logical(npoints)
    maskloop[ idxloop ] = TRUE

    generator[ maskloop ]   = 0     # near 0 to exactly 0

    maskneg = (generator < 0)   # & ! maskloop
    maskpos = (0 < generator)   # & ! maskloop

    mess    = sprintf( "%d -- %d negative, %d positive, and %d loops.\n",
                        npoints, sum(maskneg), sum(maskpos), sum(maskloop) )
    cat( "generators:       ", mess )

    #tmin    = sum( generator[maskneg] )
    #tmax    = sum( generator[maskpos] )


    cat( '\n' )
    df  = dfminmax( generator )
    mess    = sprintf( "[%g,%g]\n", df$value[1], df$value[2] )
    cat( "segment:                  ", mess )
    print( df )

    cat( '\n' )
    df  = dftrans2( generator )
    mess    = sprintf( "[%g,%g]\n", df$value[1], df$value[2] )
    cat( "2-transition subsegment:  ", mess )
    print( df )

    cat( '\n' )
    cat( "matroid:\n" )
    print( x$matroid )
    
    return( invisible(TRUE) )
    }


#   x   a zonoseg object
#   z   a M-vector of real numbers, all should be inside the zonoseg segment
#
#   returns a data frame with columns:
#       z       the given vector z
#       pcube   N x M matrix whose rows are points in the unit cube, that map to the given numbers
#               N is the # of generators 
#
#   this version tries to be as economical as possible, with the fewest number of non-zero coefficients
#   the actual code is more verbose, and not "slick"

invert.zonoseg <- function( x, z, tol=0, ... )
    {
    numpoints   = length(z)

    ok  = is.numeric(z)  &&  0<numpoints
    if( ! ok )
        {
        log_level( ERROR, "z or numpoints is invalid." )
        return(NULL)
        }

    generator   = as.double( x$matroid$matrix )
    #dfmm        = dfminmax( generator )
    #zmin        = dfmm$value[1]
    #zmax        = dfmm$value[2]

    zmin    = x$segment[1]
    zmax    = x$segment[2]

    #absgen  = abs(generator)

    idxneg  = which( generator < 0 )
    numneg  = length(idxneg)
    if( 0 < numneg )
        {
        bkneg   = c( cumsum( generator[ idxneg[numneg:1] ] )[ numneg:1 ], 0 )   #; cat( "bkneg=", bkneg, '\n' )
        }

    idxpos  = which( 0 < generator )
    numpos  = length(idxpos)
    if( 0 < numpos )    bkpos   = c( 0, cumsum( generator[idxpos] ) )


    n   = length( generator )

    pcube  = matrix( NA_real_, numpoints, n )
    colnames(pcube)    = as.character( x$matroid$ground )

    s   = numeric( n )

    for( i in 1:numpoints )
        {
        zi  = z[i]

        inside  = zmin-tol <= zi  &&  zi <= zmax+tol
        if( ! inside )
            {
            log_level( WARN, "invert.zonoseg().  zi=%g  and segment is [%g,%g].  delta=%g,%g", 
                        zi, zmin, zmax, zi-zmin, zmax-zi )
            next
            }


        if( zi == 0 )       { pcube[i, ] = 0 ; next }

        if( zi <= zmin )    { pcube[i, ] = 0 ; pcube[i,idxneg]=1 ; next }      # dfmm$pcube[1, ] ;

        if( zmax <= zi )    { pcube[i, ] = 0 ; pcube[i,idxpos]=1 ; next }       # dfmm$pcube[2, ]  }

        s[ 1:n ] = 0

        if( 0 < zi )
            {
            #   numpos must be nonzero
            j = findInterval( zi, bkpos, rightmost.closed=TRUE )

            if( 1 < j ) { s[ idxpos[1:(j-1)] ]  = 1 }

            s[ idxpos[j] ]  = (zi - bkpos[j]) / generator[ idxpos[j] ]
            }
        else if( zi < 0 )
            {
            #   numneg must be nonzero
            j = findInterval( zi, bkneg, rightmost.closed=TRUE, left.open=TRUE ) #; cat( "j=", j, '\n' )

            if( j < length(bkneg) ) { s[ idxneg[ (j+1):length(bkneg) ] ] = 1 }

            s[ idxneg[j] ]  = 1 - (zi - bkneg[j]) / (-generator[ idxneg[j] ])
            }
            
        masklo  = (s < 0)
        maskhi  = (1 < s)
        if( any(masklo | maskhi) )
            {
            count   = sum(masklo)
            if( 0 < count )
                {
                deltamax    = max( -s )
                log_level( WARN, "invert.zonoseg(). %d output values slightly less than 0 (deltamax=%g) changed to 0.",
                                        count, deltamax )
                s[ masklo ] = 0
                }
            
            count   = sum(maskhi)
            if( 0 < count )
                {
                deltamax    = max( s-1 )
                log_level( WARN, "invert.zonoseg(). %d output values slightly greater than 1 (deltamax=%g) changed to 1.",
                                        count, deltamax )
                s[ maskhi ] = 1
                }
            }

        pcube[i, ] = s
        }

    if( FALSE )
        {
        #   test it
        cat( "invert.zonoseg() test:\n" )
        print( x$matroid$matrix )
        cat( sprintf( "z=%g   zmin=%g   zmax=%g   sum(matrix)=%g\n", z, zmin, zmax, sum(x$matroid$matrix) ) )
        test    = as.double( x$matroid$matrix %*% t(pcube) )  -  z
        print( range(test) )
        }

    #rnames  = names(z)
    #if( is.null(rnames) )   rnames = 1:numpoints
    #out = data.frame( row.names=rnames )

    out = data.frame( z=z )    
    # out$z       = z
    out$pcube   = pcube

    return( out )
    }


inside_zonoseg  <- function( x, p )
    {
    ok  = is.numeric(p)  &&  0<length(p)
    if( ! ok )
        {
        log_level( ERROR, "p is invalid." )
        return(NULL)
        }

#    rnames  = names(p)
#    if( is.null(rnames)  ||  anyDuplicated(rnames)!=0 )
#        rnames  = 1:length(p)
        
    #   compute distances from x$segment endpoints
    distance    = pmax( x$segment[1] - p, p - x$segment[2] )

    out = data.frame( p=p )

    #   out$p           = p
    out$inside      = distance <= 0
    out$distance    = distance
    out$idxhyper    = 1L    # there is only 1 hyperplane, the empty set {}
    
    return( out )
    }






getsegment.zonoseg <- function( x )
    {
    return( x$segment )
    }

getsegment2trans.zonoseg <- function( x )
    {
    return( x$segment2trans )
    }

dfminmax    <- function( generator )
    {
    pcube  = rbind( ifelse(generator<0,1,0), ifelse(0<generator,1,0) )
    #colnames(pcube)   = names(generator)

    #pcube  = rbind( pcube, 0.5 )

    out = data.frame( row.names=c('zmin','zmax') )
    out$value   = pcube %*% generator
    out$pcube  = pcube

    return( out )
    }

dftrans2    <- function( generator )
    {
    # treat special cases
    for( pass in 1:2 )
        {
        if( pass == 1 )
            mask    = (0 < generator)
        else
            mask    = (0 <= generator)  # only different if generator has a 0

        dfrun   = findRunsTRUE(mask,TRUE)
        if( nrow(dfrun) <= 1 )
            {
            #   mask is optimal and has 2 (or 0) transitions, so we are done with little work
            source_max  = as.numeric(mask) #; print(source_max  )

            source_min  = 1 - source_max    #; print( source_min )

            source  = rbind( source_min, source_max )

            colnames(source) = names(generator)

            out = data.frame( row.names=c('zmin-2trans','zmax-2trans') )
            out$value   = source %*% generator
            out$source  = source

            return(out)
            }
        }

    #   now we have to do more work
    generator2  = rep( generator, 2 )

    #   find local mins and maxs in cumsum(generator2)
    minlocal    = integer(0)
    maxlocal    = integer(0)

    genprev = generator2[1]

    if( genprev == 0 )
        {
        #   wrap around to find a non-zero one
        idx     = which( sign(generator2) != 0 )
        genprev = generator2[ idx[ length(idx) ] ]
        }

    n   = length(generator2)

    for( k in 2:n )
        {
        gen = generator2[k]

        if( gen == 0 )    next    # ignore it

        if( genprev<0  &&  0<gen  )
            #   local min
            minlocal    = c( minlocal, k-1 )
        else if( 0<genprev  &&  gen<0  )
            #   local max
            maxlocal    = c( maxlocal, k-1 )

        genprev   = gen
        }

    #cat( "minlocal =", minlocal, '\n' )
    #cat( "maxlocal =", maxlocal, '\n' )

    minmax  = expand.grid( minlocal, maxlocal )
    colnames(minmax)    = c('minlocal','maxlocal')
    #print( minmax )

    rowdiff = minmax[ ,2] - minmax[ ,1]
    minmax  = minmax[ 0<rowdiff  &  rowdiff<n, ]    #; print( minmax )

    cs  = cumsum( generator2 ) #; print(cs)

    delta   = cs[ minmax[ ,2] ] - cs[ minmax[ ,1] ]

    minmax  = cbind( minmax, delta ) #; print( minmax )

    idx = which.max(delta)[1]
    #cat( "argmin=", minmax[idx,1], "    argmax=", minmax[idx,2], "   sum=", delta[idx], '\n' )

    imin    = minmax[idx,1] + 1L    # the 1L is tricky !   missed it first time
    imax    = minmax[idx,2]

    #   now make it wrap around if necessary
    run = ( ((imin:imax)-1L) %% length(generator) ) + 1L

    source_max  = numeric( length(generator) )
    source_max[run] = 1     #; print(source_max  )

    source_min  = 1 - source_max    #; print( source_min )

    source  = rbind( source_min, source_max )

    out = data.frame( row.names=c('tmin-2trans','tmax-2trans') )
    out$value   = source %*% generator
    out$source  = source

    return( out )
    }
    
    
#   x   a matroid, generated from a matrix
#
#   returns a list of zonoseg's
#   if there is no matrix, or no multiples, then it returns NULL

makezonoseglist <- function( x )
    {
    mat     = getmatrix( x )

    if( is.null(mat) )  return(NULL)

    nummultiples    = length(x$multiple)

    if( nummultiples == 0 )  return(NULL)

    gndorig         = getground(x)

    idxfromgroundORIG   = idxfromgroundfun( gndorig )

    out = vector( nummultiples, mode='list' )

    for( i in 1:nummultiples )
        {
        cmax    = x$multiplesupp$cmax[i]

        idxcol  = idxfromgroundORIG[ x$multiple[[i]] ]

        #cat( "cmax=", cmax, "   idxcol=", idxcol, '\n' )

        mat1        = mat[ cmax, idxcol, drop=FALSE] # ;        print( mat1 )

        out[[i]]    = zonoseg( mat1, ground=gndorig[idxcol] )
        }

    return( out )
    }

    
    

if( FALSE )
{
is_pointed.zonoseg <- function( x )
    {
    return( 0 %in% x$segment )
    }

is_salient.zonoseg <- function( x )
    {
    return( 0 %in% x$segment )
    }
}


    
if( FALSE )
{
minkowskisum.zonoseg <- function( zono1, zono2, e0=0, ground=NULL, ... )
    {
    if( ! inherits(zono2,"zonoseg") )
        {
        log_level( ERROR, "2nd argument zono2 is invalid." )
        return(NULL)
        }

    #   get the 2 matrices and cbind them
    mat1    = zono1$matroid$matrix
    mat2    = zono2$matroid$matrix
    mat = cbind(mat1,mat2)

    #gnd1    = getground( zono1$matroid )
    #gnd2    = getground( zono2$matroid )
    #gnd     = c( gnd1, gnd1[length(gnd1)] + gnd2 )

    out = zonoseg( mat, e0=e0, ground=ground )

    return( out )
    }


'%+%.zonoseg'  <- function(zono1,zono2)
    {
    return( minkowskisum( zono1, zono2 ) )
    }
}





#--------       UseMethod() calls           --------------#

getsegment <- function( x )
    {
    UseMethod("getsegment")
    }

getsegment2trans <- function( x )
    {
    UseMethod("getsegment2trans")
    }







#####################   deadwood below  ###################################


#       shift       the sum of the negative generators, subtracting this from the zonoseg
#                   translates it to a "standard" segment [0,X]

#   x   a zonoseg object
#   z   a vector of real numbers, all should be inside the zonoseg segment

invert_old.zonoseg <- function( x, z, ... )
    {
    n   = length(z)

    ok  = is.numeric(z)  &&  0<n
    if( ! ok )
        {
        log_level( ERROR, "z is invalid." )
        return(NULL)
        }

    generator   = as.double( x$matroid$matrix )
    dfmm        = dfminmax( generator )
    zmin        = dfmm$value[1]
    zmax        = dfmm$value[2]

    absgen  = abs(generator)

    bkpnt   = c( 0, cumsum(absgen) )

    maskneg = (generator < 0)

    m   = length( x$matroid$ground )

    source  = matrix( NA_real_, n, m )
    colnames(source)    = as.character( x$matroid$ground )

    s   = numeric( m )

    for( i in 1:n )
        {
        zi  = z[i]

        inside  = zmin <= zi  &&  zi <= zmax
        if( ! inside )  next

        if( zi == zmin )    { source[i, ] = dfmm$source[1, ] ; next }

        if( zi == zmax )    { source[i, ] = dfmm$source[2, ] ; next }

        zprime  = zi - x$shift

        j = findInterval( zprime, bkpnt, rightmost.closed=TRUE )

        s[ 1:m ] = 0
        if( 1 < j ) { s[ 1:(j-1) ]  = 1 }

        s[ j ]  = (zprime - bkpnt[j]) / absgen[j]

        #   negative generators are special
        s[maskneg]  = 1 - s[maskneg]

        source[i, ] = s
        }

    rnames  = names(z)
    if( is.null(rnames) )   rnames = 1:n

    out = data.frame( row.names=rnames )

    out$z       = z
    out$source  = source

    return( out )
    }

