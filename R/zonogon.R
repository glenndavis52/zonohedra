#
#   zonogon is a 2-dimensional zonotope
#
#   implemented as a list with items:
#       matroid     the matroid for the zonogon, which includes the generating matrix, etc.
#       center      a 2-vector
#       facet       a data.frame with a row for each hyperplane in the simplified matroid, and in the same order
#                   and for each facet of the zonogon, and these columns:
#               center  the center of the edge of the centered zonogon
#               normal  outward-pointing unit normal
#               beta    equation of the slab is  -beta <= <x,normal> <= beta.  We always have beta > 0.
#       facet0      integer vector of facet/hyperplane indexes that contain the point 0
#       vertex      (2N)x2 matrix with vertex in each row, not centered.  N is the number of simplified generators.
#       pcube       (2N)xN matrix of raw 0-1 values, that maps to vertex[,]
#       tilingdata  data.frame with a row for each pgram tile in the standard tiling, and these columns
#               idxpair     2 indexes of the pgram
#               center      center of the pgram tile, relative to the center of the zonogon

#   NOTE:  If the zonogon is a pgram, tilingdata has only one row !

#   The tiling has N pgrams meeting 0, where N is the number of generators in the simplified matroid.
#   and there is only 1 pgram meeting the white point.



#   zonogon constructor
#
#   mat     a numeric matrix with 1 row
#   e0      threshold for a column vector to be considered 0
#   e1      threshold for codirectionality
#   ground  ground set, an integer vector in increasing order and length(ground) = ncol(x)

zonogon <- function( mat, e0=0, e1=1.e-6, ground=NULL )
    {
    ok  = is.matrix(mat)  &&  is.numeric(mat)  &&  nrow(mat)==2  &&  2<=ncol(mat)
    if( ! ok )
        {
        log_level( ERROR,  "mat is invalid." )
        return(NULL)
        }


    matroid     = matroid( mat, e0=e0, e1=e1, ground=ground )
    if( is.null(matroid) )
        return(NULL)

    out = list()

    class( out )    = c( "zonogon", "zonotope", class(out) )

    out$matroid = matroid

    #   get the simplified generators, and their ground set
    matsimple   = getsimplified(matroid)

    matgen  = getmatrix( matsimple )
    gndgen  = getground( matsimple )

    #   at this point we use the simplifying fact
    #   that the n generators in matgen and the n hyperplanes in matsimple
    #   are in the same order

    n   = ncol(matgen)  # number of generators and hyperplanes

    #   make a data.frame with all generators and their opposites
    facet   = data.frame( idx=rep(1:n,2), sign=c(rep(1L,n),rep(-1L,n)) )

    #   compute all the outward pointing edge normals
    normal  = cbind( matgen[2, ], -matgen[1, ] )

    #   now unitize
    normal  = normalizeMatrix( normal, 1L )
    facet$normal  = rbind(normal,-normal)

    #   sort normals and their opposites in cclockwise order,
    #   starting with the smallest positive angle
    theta   = atan2( facet$normal[ ,2], facet$normal[ ,1] )
    maskneg =  theta<0
    theta[ maskneg ] = theta[ maskneg ] + 2*pi
    perm    = order( theta )

    #   but only keep the first n rows
    facet = facet[ perm[1:n], ]

    #   trace n edges in cclockwise order, initialize the 1st edge
    #   matrix with the center of each edge
    center  = matrix( NA_real_, n, 2 )

    vertex  = matrix( NA_real_, n, 2 )


    #   facet0 is the vector of all facet pairs that contain 0 and the "whitepoint"
    #   the length of facet0 is one of:
    #       0   0 is in the interior of the zonogon
    #       1   0 is in the interior of an edge of the zonogon (rare)
    #       2   0 is a vertex of the zonogon
    #   This facet0 vector depends on whether we are considering the original generators,
    #   or the generators from the matrix of the simplified matroid.
    #   We compute it initially for the simplified generators, and then modify it later.
    facet0  = integer(0)

    #   outward pointing normal for the first edge
    normal  = facet$normal[1, ]

    #   compute center of first edge
    pcube   = 0.5 * sign( as.double(normal %*% matgen) )    #; cat( "j=", facet$idx[1], "  pcube=", pcube, '\n' )

    #   ensure that coordinate face$idx[1] is *exactly* zero
    #   variable j is carried into the following for() loop, and is updated during each iteration
    j   = facet$idx[1]
    pcube[ j ]   = 0

    cubemat = matrix( NA_real_, nrow=n, ncol=length(pcube) )

    for( k in 1:n )
        {
        #j   = facet$idx[k]

        if( TRUE )
            {
            #   test that pcube really maps to the center of an edge
            idx0    = which( pcube==0 )
            if( length(idx0) != 1 )
                {
                log_level( FATAL, "Internal Error. pcube does not have exactly one 0 value. k=%d.", k )
                # cat( "pcube=", pcube, '\n' )
                return(NULL)
                }

            #   delete coord idx0 (which has value 0) and test that the others are +0.5 or -0.5
            if( ! all( abs(pcube[-idx0]) == 0.5 ) )
                {
                log_level( FATAL, "Internal Error. pcube to edge is invalid. k=%d.", k )
                # cat( "pcube=", pcube, '\n' )
                return(NULL)
                }
            }

        #   map pcube to center[k, ]
        center[k, ] = matgen %*% pcube

        if( allequalskip( pcube, j ) )
            {
            #   this edge has 0 as one of the endpoints
            facet0  = c( facet0, j )
            #cat( "k=", k, "    pcube=", pcube, '\n' )
            #cat( "facet0=", facet0, '\n' )
            }

        #   advance from edge center to vertex of zonogon
        pcube[j]    = facet$sign[k] * 0.5       # all coords of pcube are now +1/2 or -1/2

        if( TRUE )
            {
            #   test that pcube really maps to a vertex of the zonogon
            if( ! all( abs(pcube) == 0.5 ) )
                {
                log_level( FATAL, "Internal Error. pcube to vertex is invalid. k=%d.", k )
                # cat( "pcube=", pcube, '\n' )
                return(NULL)
                }
            }


        #   map pcube to vertex and save it
        vertex[k, ] = matgen %*% pcube

        #   save pcube for vertex k
        cubemat[k, ]    = pcube

        if( k < n )
            {
            #   advance j to the next edge
            #   advance pcube to the center of the next edge
            j = facet$idx[k+1]
            pcube[ j ]    = 0     #;    cat( "j=", j, "  pcube=", pcube, '\n' )
            }
        }

    if( TRUE  && !( length(facet0) %in% c(0,2) ) )
        {
        log_level( FATAL, "Internal Error. length(facet0)=%d, which is not 0 or 2.", length(facet0) )
        return(NULL)
        }

    facet$center  = center

    #   put the facets in column index order, which matches the hyperplane index order
    #   this is currently assumed by section.zonogon
    facet   = facet[ order( facet$idx ), ]

    #   drop columns idx and sign since they are no longer needed
    facet$idx   = NULL
    facet$sign  = NULL

    #   calculate the n plane constants beta
    #   these plane constants are for the centered zonogon
    facet$beta    = rowSums( facet$normal * facet$center )

    #   since matgen[] has rank 2, 0 is in the interior of the _centered_ zonogon,
    #   and this implies that all n beta's are positive.
    #   verify this
    betamin = min( facet$beta  )
    if( betamin <= 0 )
        {
        log_level( FATAL, "Internal Error.  min(beta)=%g <= 0.", betamin )
        return(NULL)
        }

    out$center  = 0.5 * .rowSums( mat, nrow(mat), ncol(mat) )

    out$facet   = facet


    #   facet0 is now for the simplified generators.
    #   modify facet0[] for the original generators, by removing some of them.
    #   get the "mixed direction" generators,
    #   which are the generator minors that are non-zero
    colidx      = getmixed( matroid )   # column indexes of the simplified matrix, but pass the original matroid

    if( 0<length(facet0)  &&  0<length(colidx) )
        {
        #   keep only those hyperplanes in facet0 that contain all mixed generators, and remove all the rest
        hyper   = matsimple$hyperplane[ facet0 ]

        #   keep is a logical vector entries corresponding to those in facet0 and hyper
        keep    = .Call( C_issuperset, hyper, gndgen[colidx] )

        #print( hyper )
        #cat( "gndgen[colidx]=", gndgen[colidx], '\n' )
        #cat( "keep= ", keep, '\n' )
        #cat( "changing facet0 from ", facet0, " to ", facet0[keep], '\n' )

        facet0  = facet0[ keep ]

        if( FALSE )
            {
            if( length(colidx) == 1 )
                {
                #   colidx is a single facet;  is colidx in facet0 ?
                k   = match( colidx, facet0 )
                if( ! is.na(colidx) )
                    #   keep colidx and throw away the other one
                    facet0  = facet0[k]
                else
                    #   facet0 should be empty, for the original generators.
                    facet0  = integer(0)
                }
            else
                #   there are 2 or more "mixed direction" generators
                #   so 0 is in the interior of the zonogon, and facet0 should be empty.
                facet0  = integer(0)
            }
        }


    out$facet0  = facet0


    #   add vertices by symmetry, so we complete the *full* boundary
    vertex  = rbind( vertex, -vertex )

    #   translate to original and non-centered zonogon
    #  out$vertex = vertex + matrix( out$center, nrow(vertex), 2, byrow=TRUE )      # back to original coords, slow
    out$vertex  = .Call( C_sumMatVec, vertex, out$center, 1L )                      # this is faster


    #   add the cube coords by symmetry too, so we complete the *full* boundary
    cubemat = rbind( cubemat, -cubemat )

    #   record cubemat which has values +/- 0.5
    #   but convert to raw bytes in 0-1 format to save memory
    out$pcube       = as.raw( 0 < cubemat )
    dim(out$pcube)  = dim(cubemat)


    #   add data for the standard tiling
    tilingdata  = NULL

    if( n == 2 )
        {
        #   trivial tiling, only 1 tile
        tilingdata          = data.frame( row.names=1 )
        tilingdata$idxpair  = matrix( 1:2, nrow=1, ncol=2 )
        tilingdata$center   = matrix( 0, nrow=1, ncol=2 )
        }
    else
        {
        #   3 or more generators,  n*(n-1)/2 tiles
        zono    = liftedzonohedron( matgen, ground=gndgen )

        if( ! is.null(zono) )
            {
            matroidsimp3    = getsimplified(zono$matroid)

            if( TRUE  &&  ! all( lengths(matroidsimp3$hyperplane)==2 ) )
                {
                log_level( FATAL, "Internal Error. Some hyperplanes of lifted zonohedron are non-trivial." )
                return(NULL)
                }

            #   make lookup table from ground to column index
            idxfromground   = idxfromgroundfun( gndgen )

            pgrams  = nrow(zono$facet)  #   # of tiles = n*(n-1)/2

            idxpair     = matrix( NA_integer_, nrow=pgrams, ncol=2 )

            for( k in 1:pgrams )
                {
                gen2            = matroidsimp3$hyperplane[[k]]
                idxpair[k, ]    = idxfromground[ gen2 ]

                #   take the center of the pgram on the zonohedron facet, and drop the Z          # and the zonogon center
                #    centerpgram[k, ]    = zono$facet$center[k,1:2]    #+ out$center
                }

            centerpgram = zono$facet$center[ ,1:2]

            if( TRUE )
                {
                #   verify the correct order
                perm    = order( idxpair[ ,1], idxpair[ ,2] )
                if( ! all( perm == 1:length(perm) ) )
                    {
                    log_level( WARN, "idxpair[,] is not in the correct order !" )
                    return(NULL)
                    }
                }

            tilingdata          = data.frame( row.names=1:pgrams )
            tilingdata$idxpair  = idxpair
            tilingdata$center   = centerpgram
            }
        }

    out$tilingdata  = tilingdata

    return(out)
    }


#   n   a positive integer, so the step size on the circle is 2*pi/n
#   m   number of points to compute, starting at 1

polarzonogon    <- function( n, m=n, ground=NULL )
    {
    if( is.null(ground) )
        ground  =  1L:m
    else if( length(ground) != m )
        {
        log_level( ERROR, "ground is invalid, because the length is incorrect." )
        return(NULL)
        }

    if( n < 3 )
        {
        log_level( ERROR, "n=%d is invalid.", n )
        return(NULL)
        }

    if( m<2 || n<m )
        {
        log_level( ERROR, "m=%d is invalid.", m )
        return(NULL)
        }

    u   = (0:(m-1)) / n

    mat = t( tocircle(2*pi*u) )

    return( zonogon(mat,ground=ground) )
    }

if( FALSE )
{
#   x   a zonogon object
#   p   an Mx2 matrix, etc.
#
#   value   see inside_zonotope()

inside.zonogon <- function( x, p )
    {
    p   = prepareNxM( p, 2 )
    if( is.null(p) )    return(NULL)

    #   translate p to the centered zonogon
    gcentered   = p - matrix( x$center, nrow(p), 2, byrow=TRUE ) #; print(gcentered)

    hg  = tcrossprod( x$facet$normal, gcentered )    #; print( str(hg) )

    out = inside_zonotope( x, p, hg )

    return( out )
    }


#   x           a zonogon object
#   direction   Mx2 matrix, with the M directions in the rows, direction (0,0) is invalid
#   tol         tolerance for argmax, being in the same affine subspace

#   value:   see support_zonotope()

support.zonogon <- function( x, direction, tol=5.e-15 )
    {
    return( support_zonotope(x,direction,tol) )
    }
}



#   x           a zonogon object
#   base        a numeric vector of length 2, the basepoint of all the rays
#               base must be in the interior of x,
#               or if x is non-negative, base can also be the black or white point on the boundary(x)
#   direction   an Nx2 matrix with directions in the rows
#
#   value   a dataframe with columns
#           base        given basepoint of the ray (all the same)
#           direction   given direction of the ray
#           facetidx    of the facet (edge) where ray exits the zonogon
#           sign        +1 if beta, and -1 if -beta
#           tmax        ray parameter of intersection with facet
#           point       point of intersection with the intersection with facet
#           timetrace   computation time (sec)

raytrace.zonogon <- function( x, base, direction, plot=FALSE, ... )
    {
    ok  = is.numeric(base)  &&  length(base)==2  &&  all( is.finite(base) )
    if( ! ok )
        {
        log_level( ERROR, "base is invalid. It must be a numeric vector of length 2, and all entries finite." )
        return(NULL)
        }

    direction   = prepareNxM( direction, 2 )
    if( is.null(direction) )    return(NULL)

    #   translate base to the centered zonogon
    base    = as.numeric(base)

    gcentered   = base - x$center  #; print(gcentered)

    dim(base)       = c(1,2)
    dim(gcentered)  = c(1,2)

    hg  = tcrossprod( x$facet$normal, gcentered )    #; print( str(hg) )

    # hg  = as.numeric( x$facet$normal  %*%  gcentered )      #; print( str(hg) )

    #df  = inside_zonotope( x, base, hg )

    #if( ! df$inside )
    #    {
    #    log_level( ERROR, "point base=(%g,%g) is not in the interior of the zonogon.",
    #                            base[1], base[2] )
    #    return(NULL)
    #    }


    dim(gcentered)  = NULL

    #   test whether base is black or white point, no tolerance here
    blackpt = ifelse( is_salient(x),  all(gcentered == -x$center), FALSE )
    whitept = ifelse( is_salient(x),  all(gcentered ==  x$center), FALSE )

    if( blackpt || whitept )
        {
        #   get the normals for all facets that meet 0
        # normal0 is Mx2 where M is the number of these facets
        normal0 = x$facet$normal[ x$facet0,  , drop=FALSE ]

        #   get a vector that points from base into the interior
        if( blackpt )
            interiorvec = x$center
        else
            interiorvec = -(x$center)

        test    = normal0 %*% interiorvec
        dim(test)   = NULL

        #   in the next line, sign(test) is replicated to all columns of normal0
        #   normal0 will be used in the for() loop below
        normal0 = sign(test) * normal0

        if( FALSE )
            {
            cat( "before test=", test, "\n" )
            test    = normal0 %*% interiorvec
            cat( "after test=", test, "\n" )
            }
        }
    else
        {
        #   not blackpt or whitept, so verify that base is in the *interior* of x
        df  = inside_zonotope( x, base, hg )

        if( 0 <= df$distance )
            {
            log_level( ERROR, "point base=(%g,%g) is not in the interior of the zonogon. distance=%g >= 0",
                                    base[1], base[2], df$distance )
            return(NULL)
            }
        }


    dim(base)       = NULL

    n   = nrow(direction)

    tmax        = rep(NA_real_,n)
    idx         = rep(NA_integer_,n)
    sign        = rep(NA_integer_,n)
    point       = matrix(NA_real_,n,2)
    timetrace   = rep(NA_real_,n)

    for( k in 1:n )
        {
        time_start  = gettime()

        v   = direction[k, ]

        if( any( is.na(v) ) )   next

        if( all( v==0 ) ) next    # 0-vector

        if( blackpt || whitept )
            {
            interior    = all( 0 < normal0 %*% v )

            if( ! interior )
                #   the ray starts on the boundary and does *not* enter the interior, so give up
                next
            }

        hv      = x$facet$normal  %*%  v

        numerator   = x$facet$beta - sign(hv)*hg
        tvec    = numerator / abs(hv)

        tvec[ ! is.finite(tvec) ]   = Inf

        j   = which.min( tvec )     # this ignores ties and Infs

        if( tvec[j] <= 0 )  next    # failed to intersect properly

        tmax[k] = tvec[j]       # tmax[k] is not negative

        idx[k]  = j  # x$facet$idx[j]

        sign[k] = as.integer( sign(hv[j]) )

        optcentered     = gcentered  +  tmax[k] * v     #; print( optcentered )
        point[k, ]      = optcentered + x$center

        timetrace[k]   = gettime() - time_start
        }

    rnames  = rownames(direction)
    if( is.null(rnames)  ||  anyDuplicated(rnames) )   rnames = 1:n

    out = data.frame( row.names=rnames )

    out$base        = matrix( base, n, 2, byrow=TRUE )  # replicate base to all rows
    out$direction   = direction
    out$facetidx    = idx
    out$sign        = sign
    out$tmax        = tmax
    out$point       = point
    out$timetrace   = timetrace

    cnames  = colnames(base)
    if( is.null(cnames) )   cnames = colnames(direction)

    colnames(out$point)   = cnames

    if( plot )
        {
        if( grDevices::dev.cur() == 1 )
            {
            log_level( WARN, "Cannot add rays to plot, because there is no plotting window open." )
            }
        else
            {
            graphics::segments( out$base[ ,1], out$base[ ,2],   point[ ,1], point[ ,2],  col='red', lty=1 )
            graphics::points( point[ ,1], point[ ,2],  col='red', pch=20 )
            }
        }


    return( out )
    }




#   section() compute intersection of zonogon and line(s)
#
#   x           a zonogon object
#   normal      a non-zero numeric vector of length 2, the normal of all the lines
#   beta        a vector of line-constants.  The equation of line k is: <x,normal> = beta[k]
#
#   value   a data.frame with length(beta) rows.  And these columns:
#           normal      the given normal and the same in each row
#           beta        given line constant of the line
#           boundary1   the 1st point, or NA
#           boundary2   the 2nd point, or NA
#
#   If there are 2 boundary points, then boundary1 is on the left and boundary2 on the right,
#   assuming that normal is "up".
#   Another way:  if normal is "north" then b1 is on "west" and b2 on the "east".

section.zonogon <- function( x, normal, beta, tol=1.e-10, plot=FALSE, ... )
    {
    ok  = is.numeric(normal)  &&  length(normal)==2  &&  all( is.finite(normal) )  &&  any( normal!=0 )
    if( ! ok )
        {
        log_level( ERROR, "normal is invalid. It must be a non-zero numeric vector of length 2, and all entries finite." )
        return(NULL)
        }

    ok  = is.numeric(beta)  &&  0<length(beta)   &&  all( is.finite(beta) )
    if( ! ok )
        {
        log_level( ERROR, "beta is invalid. It must be a numeric vector of positive length, and all entries finite." )
        return(NULL)
        }

    normal  = as.numeric(normal)

    W   = t( getmatrix( getsimplified(x$matroid) ) )

    #   compute functional in R^n
    functional  = as.numeric( W %*% normal )   #; print( functional )

    #   find vertex where functional is maximized
    vertex_max  = 0.5 * sign( functional )      #; print( vertex_max )  # this a vertex of the n-cube, translated by -1/2

    betamax = sum( functional * vertex_max )    # the maximum of <x,normal> for x in the zonotope
    betamin = -betamax  # by symmetry

    #   the facet centers are much easier to work with when the antipodal ones are added to the originals
    center      = rbind( x$facet$center, -x$facet$center )
    W2          = rbind( W, W )
    functional  = rep( functional, 2 )
    delta       = 0.5 * abs(functional)


    cn  = center %*% normal

    #   compute range of normal over each edge/facet
    cnpos   = cn + delta
    cnneg   = cn - delta

    #   translate beta to the centered zonogon
    betacentered   = as.numeric(beta) - sum( x$center * normal )  #; print(betacentered)

    n   = length(beta)

    #out         = vector( n, mode='list' )
    #names(out)  = sprintf( "normal=%g,%g. beta=%g", normal[1], normal[2],  beta )

    rnames  = names(beta)
    if( is.null(rnames)  ||  anyDuplicated(rnames)!=0 ) rnames = 1:n

    out = data.frame( row.names=rnames )
    out$normal  = matrix( normal, nrow=n, ncol=2, byrow=TRUE )
    out$beta    = as.numeric(beta)

    boundary1   = matrix( NA_real_, nrow=n, ncol=2 )
    boundary2   = matrix( NA_real_, nrow=n, ncol=2 )




    for( k in 1:n )
        {
        beta_k  = betacentered[k]
        if( beta_k < betamin-tol  ||  betamax+tol < beta_k )
            {
            # line does not intersect the zonogon. there is no section
            next
            }

        if( abs(beta_k - betamax) < tol )
            {
            #   special case - only one point of intersection
            boundary1[k, ]  = vertex_max %*% W  + x$center
            next
            }
        if( abs(beta_k - betamin) < tol )
            {
            #   special case - only one point of intersection
            boundary1[k, ]  = -vertex_max %*% W  + x$center
            next
            }

        idx = which( cnneg <= beta_k  &  beta_k < cnpos )  #; print( idx )

        count   = length(idx)
        ok  = count %in% 1:2
        if( ! ok )
            {
            # should not happen
            log_level( WARN, "count=%d unexpected.  count should be 1 or 2.", count )
            next
            }

        section = matrix( NA_real_, nrow=count, ncol=2 )

        for( i in 1:count )
            {
            j           = idx[i]
            alpha       = (beta_k - cn[j])/ functional[j]   #; print(alpha)
            section[i, ] = center[j, ] + alpha * W2[j, ] + x$center  # translate from centered back to original
            }

        if( count == 1 )
            boundary1[k, ]  = section[1, ]
        else
            {
            p1  = section[1, ]
            p2  = section[2, ]

            #   put p1 on "west" and p2 on "east"
            dif = p2 - p1   # dif must be orthogonal to normal
            test    = dif[1]*normal[2] - dif[2]*normal[1]
            if( test < 0 )
                {
                #   swap p1 and p2
                p=p1 ; p1=p2 ; p2=p
                }

            boundary1[k, ]  = p1
            boundary2[k, ]  = p2
            }


        #out[[k]]            = list()
        #out[[k]]$beta       = beta[k]
        #out[[k]]$section    = section
        #colnames(out[[k]]$section)  = names(normal)
        }

    out$boundary1   = boundary1
    out$boundary2   = boundary2

    if( plot )
        {
        if( grDevices::dev.cur() == 1 )
            {
            log_level( WARN, "Cannot add section to plot, because there is no plotting window open." )
            }
        else
            {
            graphics::segments( boundary1[ ,1], boundary1[ ,2],   boundary2[ ,1], boundary2[ ,2],  col='red', lty=2 )
            graphics::points( boundary1[ ,1], boundary1[ ,2],  col='red', pch=1 )
            graphics::points( boundary2[ ,1], boundary2[ ,2],  col='red', pch=1 )
            }
        }

    return( out )
    }



#   x       a zonogon
#   z       Nx2 matrix of points inside the zonogon
#   tol     tolerance for being inside the zonogon,
#           a very small positive number allows the points to be slightly outside
#   plot    if TRUE then plot the points on an existing plot of the zonogon
#
#   returns a data.frame with a row for each point z[i, ], and these columns

#       hyperidx    index of parallelogram that contains z[i, ]
#       hyper       Nx2 matrix of the generators of the parallelogram
#       pcube       NxM matrix of points in the original cube that maps to z[i, ]
#                   where M is the number of generators of the zonogon

invert.zonogon <- function( x, z, tol=0, plot=FALSE, ... )
    {
    z   = prepareNxM( z, 2 )
    if( is.null(z) )    return(NULL)

    df  = inside( x, z )
    if( is.null(df) )   return(NULL)

    #   print( df )


    inside  = df$distance < tol
    if( any( !inside ) )
        {
        log_level( WARN, "%d of %d points are not inside the zonogon.  Try increasing tol=%g.",
                                sum( !inside ), length(inside), tol )
        }

    n   = nrow(z)


    matroidsimp = getsimplified(x$matroid)

    matsimp = getmatrix( matroidsimp )
    gndsimp = getground( matroidsimp )

    centerdiff  = 0.5*rowSums(matsimp) - x$center

    numgen  = ncol(matsimp)

    if( 2 < numgen )
        {
        #   the usual nontrivial case
        zonoh   = liftedzonohedron( matsimp, ground=gndsimp )
        if( is.null(zonoh) )    return(FALSE)

        matroidsimp3    = getsimplified(zonoh$matroid)
        matrixsimp3     = getmatrix(matroidsimp3)

        idxfromground   = idxfromgroundfun( gndsimp )


        #   translate z to the centered zonogon
        pcentered   = .Call( C_sumMatVec, z, -x$center, 1L )



        hyperidx    = rep( NA_integer_, n )
        hyper       = matrix( NA_integer_, nrow=n, ncol=2 )

        pcube           = matrix( NA_real_, nrow=n, ncol=numgen )
        colnames(pcube) = gndsimp

        for( k in 1:n )
            {
            if( ! inside[k] )   next

            dot2    = as.double( zonoh$facet$normal[ ,1:2] %*% pcentered[k, ] )
            znormal = zonoh$facet$normal[ ,3]

            tmin    = (zonoh$facet$beta - dot2) / znormal
            imin    = which.max( tmin )
            zmin    = tmin[imin]

            #   tmax    = (-zonoh$facet$beta - dot2) / znormal
            #   zmax    = min(tmax)
            #   if( zmax < zmin )   next    # z is not inside the zonogon

            hyperidx[k] = imin
            hyper[k, ]  = matroidsimp3$hyperplane[[imin]]

            facetnormal = zonoh$facet$normal[imin, ]

            pc      = (sign( facetnormal %*% matrixsimp3 ) + 1) / 2
            dim(pc) = NULL

            #   force the generators of this parallelogram to 0
            colidx  = idxfromground[ hyper[k, ] ]
            pc[colidx]  = 0

            #   now map pc to the zonogon, giving the origin of the parallelogram
            org         = as.double(matsimp %*% pc) - centerdiff    #; cat( "org=", org, '\n' )
            mat2x2      = matsimp[ , colidx ]
            pc[colidx]  = as.double( solve( mat2x2, z[k, ] - org) )

            pcube[k, ]  = pc
            }
        }
    else
        {
        #     numgen == 2 is a trivial special case

        #   translate z to the simplified parallelogram
        psimple = .Call( C_sumMatVec, z, centerdiff, 1L )

        pcube   = t( solve( matsimp, t(psimple) ) )

        #cat( "psimple=", as.double(psimple), '\n' )
        #cat( "centerdiff=", centerdiff, '\n' )

        hyperidx    = 1L
        hyper       = matrix( gndsimp, nrow=nrow(pcube), ncol=2, byrow=T )
        }

    #   clamp pcube to [0,1], to eliminate very slightly outside points, such as -1.e-12
    pcube   = pmin( pmax( pcube, 0), 1 )


    if( plot )
        {
        if( grDevices::dev.cur() == 1 )
            {
            log_level( WARN, "Cannot add points to plot, because there is no plotting window open." )
            }
        else
            {
            points( z[ ,1], z[ ,2], col='red', pch=4 )
            }
        }

    #   lift to the original cube
    #   since pcube[,] has already been clamped, we can leave tol=0
    pcube   = invertcubepoints( x, pcube )
    if( is.null(pcube) )    return(NULL)

    if( FALSE )
        {
        #   test the inversion
        matorig = getmatrix( x$matroid )

        delta   = pcube %*% t(matorig)  -  z
        #cat( "range(delta)=", range(delta), '\n' )
        #print( t(matorig) )
        #print( delta )

        if( any( 5.e-15 < abs(delta) ) )
            {
            log_level( WARN, "Inversion test failed.  max(delta)=%g", max(abs(delta)) )
            }
        }

    rnames  = rownames(z)
    if( is.null(rnames)  ||  anyDuplicated(rnames)!=0 )   rnames = 1:n

    out = data.frame( row.names=rnames )
    out$z           = z
    out$pcube       = pcube
    out$hyper       = hyper
    out$hyperidx    = hyperidx

    return( out )
    }


print.zonogon  <-  function( x, ... )
    {
    cat( "zonogon:\n" )

    pairs   = nrow(x$facet)
    cat( "number of facets:                 ", 2*pairs, '  [', pairs, "antipodal facet-pairs]\n" )

    #pairs   = length( x$matroid$multiple )
    #cat( "facets with multiple generators:  ", 2*pairs, '\n' )

    gndsimp = getground( getsimplified(x$matroid) )

    idx     = getmixed(x$matroid)       #;which( 0 < rowSums( abs(x$matroid$multiplesupp$minor) ) )
    gens    = length(idx)
    cat( "generators with mixed-directions: ", gens, '   {', gndsimp[idx], '}\n' )

    cat( "center:                           ", x$center, '\n' )

    cat( "facets that contain 0:            (", length(x$facet0), " facets) {", paste(x$facet0,collapse=' '), "}\n" )  #, sep=''
    cat( "pointed:                          ", is_pointed(x), '\n' )
    cat( "salient:                          ", is_salient(x), '\n' )

    res = getmetrics(x)
    cat( "perimeter:                        ", res$perim, '\n' )
    cat( "area:                             ", res$area, '\n' )


    cat( '\n' )
    cat( "matroid:\n" )
    print( x$matroid )

    return( invisible(TRUE) )
    }


plot.zonogon  <- function( x, orientation=TRUE, normals=FALSE, elabels=FALSE,
                                tiling=FALSE, tlabels=FALSE,
                                trans2=FALSE, trans2type='both', ... )
    {
    xlim    = range( x$vertex[ ,1] )
    ylim    = range( x$vertex[ ,2] )

    if( normals )
        {
        xlim    = xlim + c(-1,1)
        ylim    = ylim + c(-1,1)
        }
    else if( elabels )
        {
        xlim    = xlim + 0.1*c(-1,1)
        ylim    = ylim + 0.1*c(-1,1)
        }


    plot( xlim, ylim, type='n', las=1, asp=1, lab=c(10,10,7), xlab='x', ylab='y' )
    grid( lty=1 )
    abline( h=0, v=0 )
    col     = ifelse( tiling, NA, 'white' )
    border  = ifelse( orientation ||  tiling, NA, 'black' )
    polygon( x$vertex[ ,1], x$vertex[ ,2], col=col, border=border )


    n   = nrow(x$facet)

    matroidsimp = getsimplified(x$matroid)

    matsimp = getmatrix( matroidsimp )
    gndsimp = getground( matroidsimp )
    numgen  = length(gndsimp)


    if( elabels )
        {
        cex = 0.6
        h   = strheight( "1", cex=cex )

        head    = x$facet$center  +  h * x$facet$normal
        head    = .Call( C_sumMatVec, rbind( head, -head ), x$center, 1L )

        text( head[ ,1], head[ ,2], as.character(rep(gndsimp,2)), cex=cex, col='black' )
        }


    if( tiling  &&  ! is.null(x$tilingdata) )
        {
        #   draw each parallelogram as polygon

        #   quadmat = matrix( 0, nrow=4, ncol=2 )

        edgecoeff   = matrix( c( -0.5,-0.5, -0.5,0.5, 0.5,0.5, 0.5,-0.5), 2, 4 )    # 2x4

        for( k in 1:nrow(x$tilingdata) )
            {
            col2    = x$tilingdata$idxpair[k, ]

            edge    = matsimp[ , col2 ]     # 2x2

            #quadmat[1, ]    =  -0.5 * edge[ , 1] - 0.5*edge[ , 2]
            #quadmat[2, ]    =  -0.5 * edge[ , 1] + 0.5*edge[ , 2]
            #quadmat[3, ]    =   0.5 * edge[ , 1] + 0.5*edge[ , 2]
            #quadmat[4, ]    =   0.5 * edge[ , 1] - 0.5*edge[ , 2]

            quadmat = edge %*% edgecoeff

            #   add the center of the pgram relative to the zonogon center, and the zonogon center
            centerpgram = x$tilingdata$center[k, ]  +  x$center

            polygon( quadmat[1, ] + centerpgram[1], quadmat[2, ] + centerpgram[2], border='red' )

            if( tlabels )
                {
                gen2    = gndsimp[ col2 ]
                lab = sprintf( "%d,%d", gen2[1], gen2[2] )
                text( centerpgram[1], centerpgram[2], lab, col='red', cex=0.5 )
                }
            }
        }



    if( is.logical(trans2) && trans2 )  trans2 = c(0L,n)

    if( length(trans2) == 2 )
        {
        subcomplex  = trans2subcomplex( n, trans2, trans2type )

        if( is.null(subcomplex) )   return(FALSE)

        vertexcube  = vertexfromcode( n, subcomplex$vertex[ ,1], subcomplex$vertex[ ,2] )

        vertexplane = vertexcube %*% t(matsimp)

        #centermat   = matrix( x$center, nrow(vertexplane), ncol=2, byrow=TRUE )

        xy1 = vertexplane[ subcomplex$edge[ ,1], ]      #+ centermat
        xy2 = vertexplane[ subcomplex$edge[ ,2], ]      #+ centermat

        thecolor    = 'blue'

        segments( xy1[ ,1]+x$center[1], xy1[ ,2]+x$center[2],  xy2[ ,1]+x$center[1], xy2[ ,2]+x$center[2], col=thecolor )

        wp  = 2 * x$center

        points( c(0,wp[1]), c(0,wp[2]), pch=20, cex=1.5, col=thecolor )
        }




    if( orientation )
        {
        #   draw arrows showing orientation of the generators
        matgen  = getmatrix( x$matroid )  #   generators in the columns

        if( is_simple(x$matroid) )
            {
            for( k in 1:n )
                {
                #idx     = x$facet$idx[k]         # raw column index
                # idxgnd  = x$facet$idxgnd[k]      # ground set index
                gen     = matgen[ ,k]

                center  = x$facet$center[k, ]

                for( s in c(1,-1) )
                    {
                    p0  = s*center - gen/2 + x$center
                    p1  = s*center + gen/2 + x$center
                    arrows( p0[1], p0[2], p1[1], p1[2], angle=15 )
                    points( p0[1], p0[2], pch=20 )
                    }
                }
            }
        else
            {
            #   make lookup table from ground to column index
            #idxfromground   = integer( max(x$matroid$ground) )
            #idxfromground[ x$matroid$ground ] = 1:ncol(matgen)
            #ground  = getground( getsimplified(x$matroid) )

            for( k in 1:n )
                {
                idxgnd  = gndsimp[k]     # ground set integer

                idxgrp  = getmultipleindex( x$matroid, idxgnd )     #; print( idxgrp )

                dominant    = matsimp[ , k ]

                if( 0 < idxgrp )
                    {
                    minor   = x$matroid$multiplesupp$minor[ idxgrp, ]

                    #   res = dominantDirection( matgen[ , idxfromground[idxgrp], drop=FALSE ] )
                    #print( res )

                    #   in this quotient, any vector norm will do, so use L^inf
                    beta    = max(abs(minor)) / max(abs(dominant))
                    }
                else
                    beta    = 0


                #   code==2 then arrowhead on p1 end
                #   code==3 then arrowhead on both ends
                code    = ifelse( 0<beta && beta<1, 3, 2 )

                center  = x$facet$center[k, ]

                for( s in c(1,-1) )
                    {
                    p0  = s*center - dominant/2 + x$center
                    p1  = s*center + dominant/2 + x$center

                    arrows( p0[1], p0[2], p1[1], p1[2], angle=15, code=code )

                    #   a point marking origin of both major and minor
                    q   = (1-beta)*p0 + beta*p1
                    points( q[1], q[2], pch=20 )
                    }
                }
            }
        }

    if( normals )
        {
        tail    = x$facet$center
        head    = x$facet$center  +  x$facet$normal

        tail    = .Call( C_sumMatVec, rbind( tail, -tail ), x$center, 1L )
        head    = .Call( C_sumMatVec, rbind( head, -head ), x$center, 1L )

        arrows( tail[ ,1], tail[ ,2], head[ ,1], head[ ,2], angle=15 )

        points( tail[ ,1], tail[ ,2], cex=.5, pch=21, bg='white' )
        }


    points(  x$vertex[ ,1], x$vertex[ ,2], col='black', pch=20 )

    #   plot the center
    points( x$center[1], x$center[2], cex=1, pch=21, bg='white' )

    main    = sprintf( "zonogon with %d generators and %d facets\n center=(%g,%g)",
                            length(getground(x$matroid)), 2L*n,
                            x$center[1], x$center[2] )

    if( ! is.null( attr(x,"sphering") ) )
        main    = paste( main, "  [spherized]" )

    title( main=main, cex.main=1 )

    return( invisible(TRUE) )
    }




if( FALSE )
{
#   methods taken from:
#       Optimal Whitening and Decorrelation
#       Agnan Kessy1, Alex Lewin, and Korbinian Strimmer (2016)
#
#   returns a new zonogon, as spherical as possible

spherize.zonogon <- function( x, method="ZCA", ... )
    {
    return( spherize_zonotope( x, method=method ) )
    }



#   x   a zonogon whose matroid is simple
#
#   replace each generator g by the pair -g/2 , +g/2

symmetrize.zonogon <- function( x )
    {
    if( ! is_simple(x$matroid) )
        {
        log_level( ERROR, "matroid is not simple" )
        return(NULL)
        }

    matgen  = getmatrix(x$matroid)

    matgen = cbind( -matgen/2, matgen/2 )

    gndgen  = getground(x$matroid)

    gndgen  = c( gndgen, gndgen[ length(gndgen) ] + gndgen )

    out = zonogon( matgen, e0=0, e1=0, ground=gndgen )

    return(out)
    }
}



#   x   a zonogon
#   W   a 2x2 invertible matrix.  Put matrix on the left, and vector on the right.

lintransform.zonogon <- function( x, W )
    {
    if( length(W) == 1 )
        W = W * diag(2)

    ok  = all( dim(W) == c(2,2) )
    if( ! ok )
        {
        log_level( ERROR, "argument W is invalid." )
        return(NULL)
        }

    #   verify W is OK
    res = try( solve(W), silent=TRUE )
    if( class(res)[1] == "try-error" )
        {
        log_level( ERROR, "matrix W is not invertible." )
        return(NULL)
        }
    Winv    = res

    #   just copy from x to out, and then make selective changes !
    out = x

    out$matroid = lintransform( x$matroid, W )

    #out$matroid$matrix  = W %*% x$matroid$matrix
    #if( ! is.null(out$matroid$simplified) )
    #    out$matroid$simplified$matrix  = W %*% x$matroid$simplified$matrix

    out$center  =   as.double( x$center %*% t(W) )

    out$facet$center    = x$facet$center %*% t(W)
    #out$facet$gen       = x$facet$gen %*% t(W)

    #   compute all the outward pointing edge normals
    normal  = x$facet$normal %*% Winv

    #   now unitize
    normal  = normalizeMatrix( normal, 1L )
    out$facet$normal  = normal

    #   calculate the n plane constants beta
    #   these plane constants are for the centered zonogon
    out$facet$beta    = .rowSums( normal * out$facet$center, nrow(normal), ncol(normal) )

    out$vertex  = x$vertex %*% t(W)

    out$tilingdata$center   = x$tilingdata$center %*% t(W)

    attr( out, "lintransform" ) = W

    return( out )
    }

if( FALSE )
{
is_pointed.zonogon <- function( x )
    {
    return( length(x$facet0) == 2 )
    }

is_salient.zonogon <- function( x )
    {
    return( 0 < length(x$facet0) )
    }
}

getmetrics.zonogon <- function( x )
    {
    matgen  = getmatrix( getsimplified(x$matroid) )

    lengen  = sqrt( colSums( matgen^2 ) )

    out             = list()
    out$vertices    = 2L*ncol(matgen)
    out$perim       = 2*sum(lengen)
    out$area        = sum(lengen * x$facet$beta)

    return( out )
    }


###################     deadwood below      ###########################

if( FALSE )
{
minkowskisum.zonogon  <- function( zono1, zono2, e0=0, e1=1.e-6, ground=NULL, ... )
    {
    if( ! inherits(zono2,"zonogon") )
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

    out = zonogon( mat, e0=e0, e1=e1, ground=ground )

    return( out )
    }

'%+%.zonogon'  <- function(zono1,zono2)
    {
    return( minkowskisum( zono1, zono2 ) )
    }
}




if( FALSE )
    {
    #   get the simplified generators, and their ground set
    matgen  = getmatrix( getsimplified(x$matroid) )
    #gndgen  = getground( getsimplified(x$matroid) )

    #   get the column indexes of generators that "meet" blackpt or whitept
    #facet0idx   = match( x$facet0, gndgen )

    facet0idx   = x$facet0

    matsub      = matgen[ , facet0idx  ]

    if( length(facet0idx) == 2 )
        {
        #   testmat is 2x2
        testmat = solve( matsub )
        if( whitept )   testmat = -testmat
        }
    else if( length(facet0idx) == 1 )
        {
        #   testmat is 1x2.  It is the inward pointing normal on boundary.
        testmat = matrix( c( matsub[2], -matsub[1] ), nrow=1 )

        if( blackpt )
            test    = testmat %*% (-x$center)
        else
            test    = testmat %*% x$center

        if( 0 < test )  testmat = -testmat
        }

    # cat( "testmat:\n" ) ; print( testmat )
    }


if( FALSE )
    {
        zono    = liftedzonohedron( matsimp, ground=gndsimp )

        if( ! is.null(zono) )
            {
            #   draw each parallelogram as polygon
            matroidsimp3    = getsimplified(zono$matroid)
            matsimp3        = getmatrix( matroidsimp3 )

            #   make lookup table from ground to column index
            idxfromground   = integer( gndsimp[numgen] )
            idxfromground[ gndsimp ] = 1:numgen

            quadmat = matrix( 0, nrow=4, ncol=2 )

            for( k in 1:nrow(zono$facet) )
                {
                gen2    = matroidsimp3$hyperplane[[k]]
                col2    = idxfromground[ gen2 ]

                edge    = matsimp3[ 1:2, col2 ]

                quadmat[1, ]    =  -0.5 * edge[ , 1] - 0.5*edge[ , 2]
                quadmat[2, ]    =  -0.5 * edge[ , 1] + 0.5*edge[ , 2]
                quadmat[3, ]    =   0.5 * edge[ , 1] + 0.5*edge[ , 2]
                quadmat[4, ]    =   0.5 * edge[ , 1] - 0.5*edge[ , 2]

                #   add the center of the facet (drop the Z) and the zonogon center
                center  = zono$facet$center[k,1:2] + x$center

                polygon( quadmat[ ,1] + center[1], quadmat[ ,2] + center[2], border='red' )

                if( tlabels )
                    {
                    lab = sprintf( "%d,%d", gen2[1], gen2[2] )
                    text( center[1], center[2], lab, col='red', cex=0.5 )
                    }
                }
            }
    }