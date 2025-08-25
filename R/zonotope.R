
#   zonotope.R
#
#   functions common to zonoseg, zonogon, and zonohedron


getmatrix.zonotope <- function( x )
    {
    return( x$matroid$matrix )
    }

getmatroid.zonotope <- function( x )
    {
    return( x$matroid )
    }

getcenter.zonotope <- function( x )
    {
    return( x$center )
    }


minkowskisum.zonotope  <- function( zono1, zono2, e0=0, e1=1.e-6, e2=1.e-10, ground=NULL, ... )
    {
    #   get the 2 matrices and cbind them
    mat1    = zono1$matroid$matrix
    mat2    = zono2$matroid$matrix

    m1  = nrow(mat1)
    m2  = nrow(mat2)

    if( m1 != m2 )
        {
        log_level( ERROR, "dimension mismatch m1 = %d  !=  %d = m2.", m1, m2 )
        return(NULL)
        }

    mat = cbind(mat1,mat2)
    colnames(mat) = NULL

    if( is.null(ground) )
        {
        #   concat g1 and translated g2
        g1  = getground(zono1$matroid)
        g2  = getground(zono2$matroid)

        ground  = c( g1, g1[ length(g1) ] + g2 )
        }


    if( m1 == 3 )
        out = zonohedron( mat, e0=e0, e1=e1, e2=e2, ground=ground )
    else if( m1 == 2 )
        out = zonogon( mat, e0=e0, e1=e1, ground=ground )
    else if( m1 == 1 )
        out = zonoseg( mat, e0=e0, ground=ground )
    else
        {
        log_level( ERROR, "dimension m1 = %d is invalid.", m1 )
        return(NULL)
        }

    return( out )
    }


'%+%.zonotope'  <- function(zono1,zono2)
    {
    return( minkowskisum( zono1, zono2 ) )
    }


is_pointed.zonotope <- function( x )
    {
    m   = nrow( getmatrix(x) )

    if( m == 3 )
        #   zonohedron
        return(3 <= length(x$facet0))
    else if( m == 2 )
        #   zonogon
        return(length(x$facet0) == 2)
    else if( m == 1 )
        #   zonoseg
        return(0 %in% x$segment)

    log_level( ERROR, "m = %d is invalid.", m )

    return( NULL )
    }



is_salient.zonotope <- function( x )
    {
    m   = nrow( getmatrix(x) )

    if( m == 3 )
        #   zonohedron
        return( 0 < length(x$facet0) )
    else if( m == 2 )
        #   zonogon
        return( 0 < length(x$facet0) )
    else if( m == 1 )
        #   zonoseg
        return(0 %in% x$segment)

    log_level( ERROR, "m = %d is invalid.", m )

    return( NULL )
    }


#   support.zonotope()
#
#   x           a zonogon or zonohedron object
#   direction   MxD matrix, with the M directions in the rows, direction (0,0) is invalid
#   tol         tolerance for argmax, being in the same affine subspace
#
#   returns a data.frame with M rows and these columns:
#       direction   the given matrix of directions
#       value       the value of the support function of x, in the given direction
#       argmax      a point on the boundary of x where the max is taken
#       dimension   of the set argmax, 0 means a vertex and 1 means an edge, and 2 means a face (of zonohedron)

support.zonotope <- function( x, direction, tol=5.e-15 )
    {
    genmat  = getmatrix( getsimplified(x$matroid) )

    direction   = prepareNxM( direction, nrow(genmat) )
    if( is.null(direction) )    return(NULL)

    n   = ncol(genmat)

    m   = nrow(direction)

    value       = rep( NA_real_, m )
    argmax      = matrix( NA_real_, m, nrow(genmat) )
    dimension   = rep( NA_integer_, m )

    functionalmat   = direction %*% genmat  # mxn

    if( 0 < tol )
        #   make entries close to 0, exactly 0
        functionalmat[ which( abs(functionalmat) <= tol, arr.ind=TRUE ) ] = 0

    for( i in 1:m )
        {
        functional  = functionalmat[i, ]

        fun0    = (functional==0)

        if( all(fun0) )    next    # ignore the 0 functional

        pcube  = 0.5*sign(functional)       # cube here is [-1/2,1/2]^N

        z   = as.double(genmat %*% pcube) + x$center  # z is in the non-centered zonogon

        value[i]        = sum( direction[i, ] * z )
        argmax[i, ]     = z
        dimension[i]    = sum( fun0 )
        }

    rnames  = rownames(direction)
    if( is.null(rnames)  ||  anyDuplicated(rnames)!=0 ) rnames = 1:m

    out = data.frame( row.names=rnames )
    out$direction   = direction
    out$value       = value
    out$argmax      = argmax
    out$dimension   = dimension

    return(out)
    }


#   x   a zonotope whose matroid is simple
#
#   replace each generator g by the pair -g/2 , +g/2

symmetrize.zonotope <- function( x, e0=0, e1=1.e-6, e2=1.e-10, ... )
    {
    if( ! is_simple(x$matroid) )
        {
        log_level( ERROR, "matroid is not simple" )
        return(NULL)
        }

    matgen  = getmatrix(x$matroid)

    matgen = cbind( matgen/2, -matgen/2 )
    colnames(matgen) = NULL

    gndgen  = getground(x$matroid)

    gndgen  = c( gndgen, gndgen[ length(gndgen) ] + gndgen )

    m   = nrow( matgen )

    if( m == 3 )
        out = zonohedron( matgen, e0=e0, e1=e1, e2=e2, ground=gndgen )
    else if( m == 2 )
        out = zonogon( matgen, e0=e0, e1=e1, ground=gndgen )
    else if( m == 1 )
        out = zonoseg( matgen, e0=e0, ground=gndgen )
    else
        {
        log_level( ERROR, "dimension m = %d is invalid.", m )
        return(NULL)
        }

    return(out)
    }


#   methods taken from:
#       Optimal Whitening and Decorrelation
#       Agnan Kessy1, Alex Lewin, and Korbinian Strimmer (2016)
#
#   returns a new zonogon or zonohedron, as spherical as possible

spherize.zonotope <- function( x, method="ZCA", ... )
    {
    M   = x$matroid$rank

    if( M == 1 )
        {
        log_level( WARN, "Cannot spherize a zonoseg, returning the original zonoseg." )
        return( x )
        }

    full    = c("ZCA","PCA-COR")
    idx = pmatch( toupper(method), full )

    ok  = is.finite(idx)  &&  length(idx)==1
    if( ! ok )
        {
        log_level( ERROR, "method='%s' is invalid.", method )
        return(NULL)
        }
    method  = full[idx]


    center  = x$facet$center
    center  = rbind( center, -center )  # add the other half, by symmetry

    res     = base::svd( center, nu=0, nv=M )

    if( ! all( 0 < res$d ) )
        {
        log_level( ERROR, "Internal error. center matrix is invalid." )
        return(NULL)
        }

    if( method == "ZCA" )
        {
        #   calculate the "whitening", or "sphering" matrix, which is MxM
        W   = res$v %*% diag( 1/res$d ) %*% t(res$v)
        }
    else if( method == "PCA-COR" )
        {
        sigma   = t(center) %*% center
        e       = diag( sigma ) #; print(e)

        P   = diag( 1/sqrt(e) ) %*% sigma %*% diag( 1/sqrt(e) )
        #cat( "P=", P, '\n' )

        decomp  = eigen( P, symmetric=TRUE ) #;        print( decomp )

        G       = decomp$vectors
        theta   = decomp$values

        Ptest   = G %*% diag(theta) %*% t(G)
        #cat( "Ptest=", P, '\n' )

        W   = diag( sqrt(1/theta) ) %*% t(G) %*% diag( sqrt(1/e) )
        #W   = t(W)
        }

    #   test the W
    #centerp = center %*% t(W)      #; print( t(centerp) %*% centerp )

    out = lintransform( x, W )

    attr( out, "sphering" ) = W

    return( out )
    }




#--------       UseMethod() calls           --------------#

spherize <- function( x, ... )
    {
    UseMethod("spherize")
    }

symmetrize <- function( x, ... )
    {
    UseMethod("symmetrize")
    }

lintransform <- function( x, W )
    {
    UseMethod("lintransform")
    }

is_pointed <- function( x )
{
    UseMethod('is_pointed')
}

is_salient <- function( x )
{
    UseMethod('is_salient')
}

inside <- function( x, p )
    {
    UseMethod("inside")
    }

raytrace <- function( x, base, direction, ... )
    {
    UseMethod("raytrace")
    }

section <- function( x, normal, beta, ... )
    {
    UseMethod("section")
    }


invert <- function( x, z, ... )
    {
    UseMethod("invert")
    }

invertboundary <- function( x, point, tol )
    {
    UseMethod("invertboundary")
    }

support <- function( x, direction, tol=5.e-15 )
    {
    UseMethod("support")
    }

'%+%' <- function(zono1,zono2)
    {
    UseMethod('%+%')
    }

minkowskisum <- function(zono1,zono2,...)
    {
    UseMethod("minkowskisum")
    }

getnormal <- function( x, ... )
    {
    UseMethod("getnormal")
    }

getmetrics <- function( x )
    {
    UseMethod("getmetrics")
    }


getmatrix  <- function( x )
{
    UseMethod('getmatrix')
}

getmatroid <- function( x )
    {
    UseMethod("getmatroid")
    }

getcenter <- function( x )
    {
    UseMethod("getcenter")
    }




canonicalboundary <- function( x, gndpair, cube=FALSE )
    {
    UseMethod("canonicalboundary")
    }



#   x   a zonotope object
#   p   an Mx3 matrix, etc.
#
#   value   see inside_zonotope()

inside.zonotope <- function( x, p )
    {
    m   = nrow( getmatrix(x) )

    if( m == 1 )
        {
        #   a zonoseg is special
        return( inside_zonoseg(x,p) )
        }

    p   = prepareNxM( p, m )
    if( is.null(p) )    return(NULL)

    #   translate p to the centered zonohedron
    #   subract x$center from every row of p    
    # gcentered   = p - matrix( x$center, nrow(p), m, byrow=TRUE ) #; print(gcentered)
    
    #gcentered  = duplicate( p )
    #res = .Call( C_plusEqual, gcentered, -x$center, 1L )   # gcentered point_centered in place
    #if( is.null(res) )  return(NULL)
    
    gcentered   = .Call( C_sumMatVec, p, -x$center, 1L )
    
    hg  = tcrossprod( x$facet$normal, gcentered )    #; print( str(hg) )

    out = inside_zonotope( x, p, hg )

    return( out )
    }


#   inside_zonotope()
#
#   x           a zonogon or zonohedron object
#   p           Mx2 matrix, with the M query points in the rowSums
#   gcentered   Mx2 matrix the same as p, but after subtracting the center of x$center
#   hg          NxM matrix, where N is the number of facets/hyperplanes of getsimplified(x).
#               hg = x$facet$normal %*% t(gcentered) = tcrossprod( x$facet$normal, gcentered )
#                                                                   Nx2      %*%    2xM
#
#   value   a dataframe with columns
#           p           the given Mx2 input matrix
#           inside      TRUE means inside the closed zonotope - interior OR boundary
#           distance    numeric signed distance to zonogon boundary
#                       negative or 0 means inside, and positive means in the exterior
#                       NOTE: if positive then the distance is only approximate.
#           idxhyper    the index, in the simplified matroid, of the critical facet

inside_zonotope <- function( x, p, hg )
    {
    distance    = abs(hg) - matrix( x$facet$beta, nrow(hg), ncol(hg) )   #; print(distance)

    if( TRUE )
        {
        #   this is the fastest method I found
        data        = .Call( C_whichMaxMatrix, distance, 2L )
        idx         = data[[1]]     #; print( idx )
        distance    = data[[2]]     #; print( distance )
        }
    else if( FALSE )
        {
        dtrans      = t(distance)
        idx         = max.col(dtrans)
        distance    = dtrans[ cbind( 1:nrow(dtrans), idx ) ]
        }
    else
        {
        myfun   <- function(z)  { idx=which.max(z) ; return( c(idx,z[idx]) ) }

        data    = apply( distance, 2, myfun )         #function(z) {suppressWarnings( max(z,na.rm=TRUE) ) } )  #;   print(distance)

        idx         = as.integer( data[1, ] )
        distance    = data[2, ]
        }


    #   distance[ ! is.finite(distance) ]   = NA_real_

    if( is_salient(x) )
        {
        #   special override for black
        black   = apply( p, 1, function(v) { isTRUE(all(v==0)) } )  #; print(black)
        if( any(black) )
            distance[black] = 0

        #   special override for white.  Fortunately multiplication by 0.5 and 2 preserves all precision.
        white   = apply( p, 1, function(v) { isTRUE(all(v==2*x$center)) } )  #; print(white)
        if( any(white) )
            distance[white] = 0
        }

    rnames  = rownames(p)
    if( is.null(rnames)  ||  anyDuplicated(rnames)!=0 )   rnames = 1:nrow(p)

    #gndgen  = getground( getsimplified(x$matroid) )

    out = data.frame( row.names=rnames )

    out$p           = p
    out$idxhyper    = idx    
    out$distance    = distance
    out$inside      = distance <= 0

    return(out)
    }




#   x       a zonogon or zonohedron
#   pmat    P x N matrix of points in the N-cube, where N is the number of generators of the original matroid
#
#   return  P x M matrix of points in the M-cube, where M is the number of generators of the simplified matroid
#           the new M-D point maps the same point in the zonogon or zonohedron as the given N-point,
#           except for the possible offset from the original zono* and the simplified zono*
#           the offset happens when there is a multiple group with "mixed" generators

projectcubepoints <- function( x, pmat )
    {
    ok  = is.matrix(pmat)  &&  is.numeric(pmat)
    if( ! ok )
        {
        log_level( ERROR, "matrix pmat is invalid." )
        return(NULL)
        }

    matorg  = getmatrix( x$matroid )

    n   = ncol(matorg)

    if( ncol(pmat) != n )
        {
        log_level( ERROR, "matrix pmat has ncol(pmat)=%d  !=  %d.", ncol(pmat), n )
        return(NULL)
        }

    if( is_simple(x$matroid) )  return( pmat )      # no change !!


    gndorg  = getground( x$matroid )

    matsimp = getmatrix( getsimplified(x$matroid) )

    m   = ncol(matsimp)

    idxfromgroundORIG    = idxfromgroundfun( gndorg )


    nummultiples    = length(x$matroid$multiple)

    if( nrow(x$matroid$multiplesupp) != nummultiples )
        {
        log_level( ERROR, "nrow(x$matroid$multiplesupp)=%d  !=  %d.",
                            nrow(x$matroid$multiplesupp), nummultiples )
        return(NULL)
        }

    multiplesupp    = x$matroid$multiplesupp

    #if( 0 < nummultiples  &&  any( multiplesupp$mixed ) )
    #    {
    #    log_level( ERROR, "Cannot project points, because the matroid has a multiple group with mixed directions." )
    #    return(NULL)
    #    }


    if( ! is.null( x$zonoseg ) )
        {
        #   zonoseg list is already computed
        #cat( "found zonoseg[[]] list.\n" )
        zonoseg = x$zonoseg
        }
    else if( 0 < nummultiples )
        {
        #   compute the zonoseg list now
        zonoseg = makezonoseglist( x$matroid )
        }

    affectedORIG = logical(n)
    affectedORIG[ idxfromgroundORIG[ x$matroid$loop ] ]                 = TRUE
    affectedORIG[ idxfromgroundORIG[ fastunion(x$matroid$multiple) ] ]  = TRUE
    #cat( "affectedORIG=", affectedORIG, '\n' )

    affectedSIMP = logical(m)
    affectedSIMP[ multiplesupp$colidx ] = TRUE
    #cat( "affectedSIMP=", affectedSIMP, '\n' )


    points  = nrow(pmat)

    out = matrix( NA_real_, points, m )
    rownames(out)   = rownames(pmat)
    colnames(out)   = colnames(matsimp)

    for( k in 1:points )
        {
        if( ! is.finite( pmat[k,1] ) ) next    # bad point

        pcube   = pmat[k, ]

        #cat( "simple pcube=", pcube, '\n' )

        #   ok  = all( 0 <= pcube  &  pcube <= 1 )

        pout    = rep( NA_real_, m )

        #   assign the unaffected coords, which we hope are in the majority
        pout[ ! affectedSIMP ] = pcube[ ! affectedORIG ]

        #   handle the multiple groups
        for( i in seq_len(nummultiples) )
            {
            #cat( "multiple ", i, '\n' )

            #cmax    = multiplesupp$cmax[i]
            mat1    = getmatrix( zonoseg[[i]] )
            #idxcol  = idxfromgroundORIG[ x$matroid$multiple[[i]] ]

            seg     = getsegment( zonoseg[[i]] )
            zmin    = seg[1]
            zmax    = seg[2]

            #   compute point z in the zonoseg.
            #cmax    = multiplesupp$cmax[i]
            idxcol  = idxfromgroundORIG[ x$matroid$multiple[[i]] ]
            #cat( "cmax=", cmax, "   idxcol=", idxcol, '\n' )

            #mat1    = matorg[ cmax, idxcol, drop=FALSE]     # ;        print( mat1 )

            z   = as.double( mat1 %*% pcube[idxcol] )       # z is just a scalar

            #   compute lambda from z
            # lambda=0 -> minor  and  lambda=1 -> major
            # lambda  = ( z  -  multiplesupp$minor[i,cmax] )  / ( multiplesupp$major[i,cmax]  -  multiplesupp$minor[i,cmax] )
            # lambda  = z  / multiplesupp$major[i,cmax]

            #   map interval [zmin,zmax] in the zonoseg to [0,1] in cube
            pout[ multiplesupp$colidx[i] ] = (z - zmin) / (zmax - zmin)

            #   cat( "z=", z, "    x =", (z-zmin)/(zmax-zmin),  "   multiplesupp$colidx[i]=", multiplesupp$colidx[i], '\n' )
            }


        #   loops in the original matroid can be ignored, they are dropped here

        out[k, ]    = pout
        }

    return( out )
    }




#   x       a zonogon or zonohedron
#   pmat    P x M matrix of points in the M-cube, where M is the number of generators of the simplified matroid
#   tol     boundary tolerance
#
#   return  P x N matrix of points in the N-cube, where N is the number of generators of the original matroid
#           the new N-vector maps the same point in the zonogon or zonohedron as the given M-point,
#           except for the offset from the original zono* and the simplified zono*

invertcubepoints <- function( x, pmat, tol=5.e-15 )
    {
    ok  = is.matrix(pmat)  &&  is.numeric(pmat)
    if( ! ok )
        {
        log_level( ERROR, "matrix pmat is invalid." )
        return(NULL)
        }

    matsimp = getmatrix( getsimplified(x$matroid) )

    if( ncol(pmat) != ncol(matsimp) )
        {
        log_level( ERROR, "matrix pmat has ncol(pmat)=%d  !=  %d.", ncol(pmat), ncol(matsimp) )
        return(NULL)
        }

    if( is_simple(x$matroid) )  return( pmat )      # no change !!


    matorg  = getmatrix( x$matroid )
    gndorg  = getground( x$matroid )

    idxfromgroundORIG    = idxfromgroundfun( gndorg )
    #cat( "idxfromgroundORIG=", idxfromgroundORIG, '\n' )

    nummultiples    = length(x$matroid$multiple)

    if( nrow(x$matroid$multiplesupp) != nummultiples )
        {
        log_level( ERROR, "nrow(x$matroid$multiplesupp)=%d  !=  %d.",
                            nrow(x$matroid$multiplesupp), nummultiples )
        return(NULL)
        }

    multiplesupp    = x$matroid$multiplesupp

    #if( 0 < nummultiples  &&  any( multiplesupp$mixed ) )
    #    {
    #    log.string( ERROR, "Cannot invert points, because the matroid has a multiple group with mixed directions." )
    #    return(NULL)
    #    }


    if( ! is.null( x$zonoseg ) )
        {
        #   zonoseg list is already computed
        #cat( "found zonoseg[[]] list.\n" )
        zonoseg = x$zonoseg
        }
    else if( 0 < nummultiples )
        {
        #   compute the zonoseg list now
        zonoseg = makezonoseglist( x$matroid )
        }

    n   = length(gndorg)

    affectedORIG = logical(n)
    affectedORIG[ idxfromgroundORIG[ x$matroid$loop ] ]                 = TRUE
    affectedORIG[ idxfromgroundORIG[ fastunion(x$matroid$multiple) ] ]  = TRUE
    #cat( "affectedORIG=", affectedORIG, '\n' )

    affectedSIMP = logical(ncol(pmat))
    affectedSIMP[ multiplesupp$colidx ] = TRUE
    #cat( "affectedSIMP=", affectedSIMP, '\n' )


    #   make simple lookup tables
    jnext   = c( 2:n, 1L )
    jprev   = c( n, 1:(n-1L) )

    points  = nrow(pmat)
    m       = ncol(pmat)


    out = matrix( NA_real_, points, n )
    rownames(out)   = rownames(pmat)
    colnames(out)   = as.character( x$matroid$ground )

    for( k in 1:points )
        {
        if( ! is.finite( pmat[k,1] ) ) next    # bad point

        pcube   = pmat[k, ]

        #cat( "simple pcube=", pcube, '\n' )

        #   ok  = all( 0 <= pcube  &  pcube <= 1 )

        pout    = numeric(n)    # rep( NA_real_, n )

        #   assign the unaffected coords, which we hope are in the majority
        pout[ ! affectedORIG ] = pcube[ ! affectedSIMP ]

        #   handle the multiple groups
        for( i in seq_len(nummultiples) )
            {
            #cat( "multiple ", i, '\n' )

            j   = multiplesupp$colidx[i]    # j is in 1:m

            s       = pcube[ j ]
            sprev   = pcube[ (j-2) %% m + 1L ]
            snext   = pcube[ (j) %% m + 1L ]

            #   compute point s in the zonoseg.  0 -> minor and 1 -> major
            # cmax    = multiplesupp$cmax[i]
            # s   = (1-lambda)*multiplesupp$minor[i,cmax]  +  lambda*multiplesupp$major[i,cmax]
            #cat( "lambda=", lambda, "    s=", s, '\n' )

            #   lift/invert s to the original cube
            # lift    = invert( zonoseg[[i]], s, tol=tol )
            #print( lift )
            #print( idxfromgroundORIG[ x$matroid$multiple[[i]] ] )

            vec = invertval( zonoseg[[i]], s, sprev, snext, tol=tol )
            if( is.null(vec) )  return(NULL)

            pout[ idxfromgroundORIG[ x$matroid$multiple[[i]] ] ] = vec    #  lift$pcube
            
            #mess    = sprintf( "from pcube[%d]=%g, added this to pout[]: ", j, s )
            #cat( mess, vec, " to pout[]\n" )
            }

        #   handle the loops, which we set to 0 or 1
        #   but to minimize the # of transitions,
        #   we prefer 0 in some cases, and 1 in other cases
        for( i in x$matroid$loop )
            {
            j   = idxfromgroundORIG[ i ]

            pp  = pout[ jprev[j] ]
            pn  = pout[ jnext[j] ]

            if( pp==0 || pp==1 )
                # on boundary of square
                pout[j] = pp
            else if( pn==0 || pn==1 )
                # also on boundary of square
                pout[j] = pn
            else
                {
                #   in interior of square, find out whether 0s or 1s are dominant
                if( sum( pout==0 ) < sum( pout==1 ) )
                    #  1s are dominant so choose 0
                    pout[j] = 0
                else
                    #  0s are dominant so choose 1
                    pout[j] = 1
                }
            }

        out[k, ]    = pout
        }
        
    if( FALSE )
        {
        resp_simp   = matsimp %*% t(pmat)
        
        resp_org    = matorg %*% t(out)
        
        delta   = resp_org - resp_simp
        cat( "|delta| = ", sum(abs(delta)), '\n' )
        print( delta )
        }
        

    return( out )
    }


#   zono    the zonoseg,
#   s       value to invert, in [0,1]
#   sprev   previous value in some hi-dimensional cube, in [0,1]
#   snext   next value in some hi-dimensional cube, in [0,1]
#   tol     tolerance passed to invert.zonoseg
#
#   returns  a point in cube of zono, that maps to (1-s)*zmin + s*zmax

invertval   <- function( zono, s, sprev, snext, tol=5.e-15 )
    {
    #   get the number of generators
    n   = length( zono$matroid$ground )

    lambda  = lambdavec( sprev, snext ) #; print(lambda)

    zmin    = getsegment(zono)[1]   # if the generators are not mixed, zmin is 0
    zmax    = getsegment(zono)[2]

    out     = numeric( n )

    if( 0 < lambda[1] )
        {
        #   this is the increasing part
        #   convert decreasing to increasing by taking complements
        lift    = invert( zono, s*zmin + (1-s)*zmax, tol=tol )
        out     = out + lambda[1] * (1 - lift$pcube)
        }

    if( 0 < lambda[2] )
        {
        #   this is the decreasing part
        lift    = invert( zono, (1-s)*zmin + s*zmax, tol=tol )
        out     = out + lambda[2] * lift$pcube
        }
        
    if( FALSE )
        {
        #   test it
        cat( "invertval() test:\n" )
        # print( x$matroid$matrix )
        cat( "lambda = ", lambda, '\n' )
        cat( sprintf( "s=%g   sprev=%g  snext=%g\n", s, sprev, snext ) )
        ztarget = (1-s)*zmin + s*zmax
        test    = as.double( zono$matroid$matrix %*% t(out) )  -  ztarget
        print( range(test) )
        }

    return( out )
    }


#   zono    a zonotope that is salient, which means that 0 is in the boundary
#
#   returns a normal vector so that the zonotope is in the non-negative closed halfspace
#   if zono is also pointed, then the zonotope is in the positive open halfspace (except for 0 itself)

supportingnormal0   <- function( zono )
    {
    matgen  = getmatrix( getsimplified(zono$matroid) )

    m   = nrow( matgen )

    if( ! is_salient(zono) )    return( rep(NA_real_,m) )

    if( m == 1 )
        {
        #   this is a zonoseg and *not* mixed, so all generators are the same sign, or 0
        return( sign( sum(matgen) ) )
        }

    #   get the normals for all facets that meet 0
    # normal0 is Fxm where F is the number of these facets
    normal0 = zono$facet$normal[ zono$facet0,  , drop=FALSE ]

    #   unfortunately we do no know whether these normals are inward are outward
    #   use the center to orient them all inward
    interiorvec = zono$center

    test    = normal0 %*% interiorvec
    dim(test)   = NULL

    #   in the next line, sign(test) is recycled to all columns of normal0
    normal0 = sign(test) * normal0      # use recyling rule

    #   normal0 now has rows that all point inward
    #   take their sum
    out = .colSums( normal0, nrow(normal0), ncol(normal0) )

    #   unitize
    out = out / sqrt( sum(out^2) )

    if( TRUE  &&  is_pointed(zono) )
        {
        #   verify the result
        test    = out %*% matgen

        if( ! all( 0 < test ) )
            {
            log_level( FATAL, "computed normal vector is invalid.  %d of %d generators failed the test.",
                                    sum(test<=0), length(test) )
            return(NULL)
            }
        }

    return( out )
    }




#   xprev, xnext    point in the unit square, not verified
#
#   returns a pair of coefficent weights that sum to 1, suitable for a convex combination
#

lambdavec <- function( xprev, xnext )
    {
    denom   = xprev*(1-xprev)  +  xnext*(1-xnext)

    lambda  = numeric(2)

    if( denom == 0 )
        {
        #   at a corner of the square
        if( xprev==0  &&  xnext==1 )
            #  limit is 1
            lambda[1]   = 1
        else if( xprev==1  &&  xnext==0 )
            # limit is 0
            lambda[1]   = 0
        else
            #   no limit at (0,0) and (1,1), so just choose something valid
            return( c(0,1) )
        }
    else
        {
        #   not a corner of the square
        lambda[1]   = ( (1-xprev) * xnext * ((1-xnext) + xprev) ) / denom
        }

    lambda[2]   = 1 - lambda[1]

    return( lambda )
    }

