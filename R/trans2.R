



#   getmetrics2trans( x, angles=TRUE, more=TRUE, tol=5.e-12 )
#
#   x       a zonohedron
#   angles  add dihedral angles of all edges
#   more    add more pgramdf columns
#   tol     tolerance for "edge-on" facets, as viewed from the center.
#           And also for the deficient shift.
#
#   returns list with:
#
#   generators          # of generators N, in the simplified matroid
#
#   pgramdf             a data frame on the pgrams in the 2-transition surface, with N*(N-1)/2 rows,
#                       computed initially by allpgramcenters2trans(), but then possibly modified.
#                       center and beta might be reversed; cross is never reversed.
#                       It has these columns:
#       idxpair         integer matrix with 2 columns i and j with 1 <= i < j <= n
#       gndpair         idxpair converted to the ground set, another integer matrix with 2 columns
#       center          real matrix with 3 columns containing the corresponding facet center,
#                       for the centered surface.  If linkingnum is negative this is reversed.
#       cross           unitized cross product of generators, never reversed
#       beta            coefficient of the plane equation of the facet
#                       If linkingnum is negative this is reversed.
#                       If the surface is starshaped, all beta are positive.
#
#   if the 2-trans surface is starshaped, then these additional columns are added
#       hyperplaneidx   index of the hyperplane that contains a congruent copy of this facet in the zonohedron
#       centermax       center of the pgram facet, relative to the center of the zonohedron
#                       for non-trivial zonogon facets of the zonohedron, this is computed from the standard tiling of the zonogon 
#       betamax         plane constant for this pgram facet
#       deficit         betamax - beta.  When coincident, it might not be exactly 0.
#       shift           distance between the 2 facet centers - centermax and center.
#                       When coincident, it should be 0, but may be very small because of truncation.
#       deficient       equal to shift >= tol.  a logical.
#       area            area of the pgram facet

#   linkingnumber   integer linking number with respect to the center, NA if undefined
#   signcount       integer 4-vector with -,0,+ and total sign counts. The total is N*(N-1)/2
#   signsurf        sign of linking number of center with the surface, +1, -1, or 0, or NA
#   starshaped      surface is starshaped, at the center, logical and often NA
#   injective       surface is injective, logical and often NA

#   if angles==TRUE, then anglesDH is added:
#   anglesDH        data frame with dihedral angle data. N*(N-1) rows and these columns:
#       pivot   integer index of the generator where dihedral angle pivots, the pivot of the "hinge" edge
#       wing1   index of the generator forming wing 1 of the "hinge"
#       wing2   index of the generator forming wing 2 of the "hinge"
#       level   the # of 1s in the level where the edge ends
#       angle   the external dihedral angle, where positive is convex and negative is concave
#       edgemid midpoint of the edge, in the *centered* polyhedron

#   if the 2-trans surface is starshaped, then these additional items are added to the output
#
#   areadeficient   sum of the areas of the deficient pgrams.  For both halves of the surface.
#   volumedeficient sum of the deficient volume, between surface and zonohedron. For both halves of the surface and zono.


getmetrics2trans <- function( x, angles=TRUE, more=TRUE, tol=5.e-12 )
    {
    if( ! inherits(x,"zonohedron") )
        {
        log_level( ERROR, "Argument x is invalid.  It's not a zonohedron." )
        return(NULL)
        }

    pgramdf   = allpgramcenters2trans( x )
    if( is.null(pgramdf) )  return(NULL)

    matsimple   = getsimplified( x$matroid )

    #centermat   = t(pgramdf$center)
    #if( ncol(centermat) != ncol(matsimple$crossprods) )
    #    {
    #    log.string( ERROR, "ncol(centermat)=%d  !=  %d=ncol(matsimple$crossprods).",
    #                        ncol(centermat) != ncol(matsimple$crossprods) )
    #    return(NULL)
    #    }

    #dotvec  = .colSums( centermat * matsimple$crossprods, nrow(centermat), ncol(centermat) )

    #   get the linking number with respect to the center
    linkingnum  = linkingnumber3( x, pgramdf )        # , , c(0,0,0) )

    #   alternate still needs work
    #linkingnum  = linkingnumber2( x )        # , , c(0,0,0) )

    betavec     = pgramdf$beta

    countneg    = sum( betavec < -tol )
    countpos    = sum( tol < betavec )
    countzero   = length(betavec) - countneg - countpos

    if( is.finite(linkingnum)  &&  abs(linkingnum) == 1 )
        {
        #   the usual case
        if( 0<countneg  &&  0<countpos )
            #   mixed signs
            starshaped  = FALSE
        else if( countzero == 0 )
            #   all negative or all positive
            starshaped  = TRUE
        else
            #   some dot products are "zero", degenerate and not mixed
            starshaped  = NA    # logical
        }
    else
        #   if the linking number is not +1 or -1, the surface cannot be starshaped
        starshaped = FALSE


    ground  = getground( matsimple )


    if( FALSE  &&  is.finite(starshaped) &&  ! starshaped  )
        {
        if( countneg <= countpos )
            mask    = betavec < 0
        else
            mask    = 0 < betavec

        gndpair         = ground[pgramdf$idxpair]
        dim(gndpair)    = c( length(gndpair)/2, 2 )

        df  = cbind( pgramdf, gndpair, betavec )
        print( df[mask, ] )
        }


    # signsurf is the linking number of the 2-transition polyhedral surface and 0
    # it is defined whenever 0 is not in the surface, even when surface has self-intersections, map is not injective
    # but we do not really have time to determine the sign accurately,
    # so choose the dominant facet sign and it will usually be correct

    #   crossprod lookup is outward when signsurf=1 and inward when signsurf=-1
    signsurf    = sign( linkingnum )

    if( is.finite(signsurf)  &&  signsurf < 0 )
        {
        #   use antipodal pgrams instead
        #   so we can compare beta with the corresponding betamax from the zonohedron,
        #   and center with centermax from the zonohedron.
        #   we only do this when starshaped is TRUE, see below
        #   the pgram normal vector stays the same,
        #   and when starshaped, it changes from inward pointing to outward pointing
        pgramdf$center = -(pgramdf$center)
        pgramdf$beta   = -(pgramdf$beta)
        }


    if( FALSE  &&  is.finite(starshaped)  &&  starshaped )
        {
        #   verify that all beta > 0
        masknonpos = pgramdf$beta <= tol
        if( any(masknonpos) )
            {
            log_level( ERROR, "internal error.  %d of %d beta coeffs are non-positive for strictly starshaped polyhedron.
                        tol=%g.", sum(masknonpos), length(masknonpos), tol )
            return(NULL)
            }
        }


    injective   = NA    # logical
    if( is.finite(linkingnum)  &&  abs(linkingnum)!=1 )
        injective   = FALSE
    else if( is.finite(starshaped)  &&  starshaped )
        injective   = TRUE

    out = list()

    out$generators      = ncol( getmatrix(matsimple) )
    out$pgramdf         = pgramdf
    out$linkingnumber   = linkingnum
    out$signcount       = c( negative=countneg, zero=countzero, positive=countpos, total=length(betavec) )
    out$signsurf        = signsurf
    out$starshaped      = starshaped
    out$injective       = injective


    if( angles  &&  is.finite(signsurf)  &&  signsurf != 0 )
        {
        #   get all the edge dihedral angles
        res = edgeangles2trans( x, signsurf )
        if( is.null(res) )  return(NULL)

        out$anglesDH  = res
        }

    if( more  &&  is.finite(starshaped)  &&  starshaped )
        {
        #   add more columns to out$pgramdf

        # beta for the zonohedron, where normal product is maximized.  The normal is outward pointing.
        betamax = x$facet$beta[ matsimple$hyperplaneidx ]

        pgrams  = nrow(pgramdf)

        if( length(betamax) != pgrams )
            {
            log_level( ERROR, "internal error. length(betamax)=%d != %d=pgrams.", length(betamax), pgrams  )
            return(FALSE)
            }

        # since the surface is starshaped, all out$pgramdf$beta > 0

        deficit     = betamax - out$pgramdf$beta    # always non-negative

        out$pgramdf$hyperplaneidx   = matsimple$hyperplaneidx
        out$pgramdf$betamax         = betamax
        out$pgramdf$deficit         = deficit

        #   x$facet$center[ , ] is the center of the zonogon facet.
        #   The next line is only correct for the trivial parallelogram facets.
        #   For a non-trivial facet, the next line sets centermax to the center of the zonogon facet,
        #   which is not the pgram center for any tiling, in general, and incorrect.
        #   This is fixed in the for() loop below.
        centermax   = x$facet$center[ matsimple$hyperplaneidx, , drop=FALSE]

        # sign modification
        # centermax   = x$facet$sign[ matsimple$hyperplaneidx ] * centermax   # recycling rule used here
        .Call( C_timesEqual, centermax, as.double(x$facet$sign[ matsimple$hyperplaneidx ]), 2L )     # multiply in place


        #   for the non-trivial zonogon facets, centermax needs special treatment

        for( k in seq_len( length(x$zonogon) ) )
            {
            zono    = x$zonogon[[k]]

            #   subgnd has M ints, where M > 2
            subgnd  = getground(zono$matroid)

            #   the the length of idx is M*(M-1)/2. the integers are in 1 : N*(N-1)/2
            idx = translateallpairs( subgnd, ground )

            masksmall   = deficit[idx] <= tol

            if( any( masksmall ) )
                {
                #   for the maximizing pgrams, with very SMALL deficit,
                #   use the 2-transition pgrams
                #   later, the shift will be computed as 0, and deficient will be FALSE
                idxsub  = idx[ masksmall ]

                #   assign only those centers for which deficit is SMALL
                centermax[ idxsub, ]    = out$pgramdf$center[ idxsub, ]
                }

            masklarge   = ! masksmall

            if( any( masklarge ) )
                {
                #   for the maximizing pgrams, with LARGE deficit
                #   use the tiling pgrams, for the standard tiling of the zonogon facet
                center3D    = gettilecenters3D( x, k )

                #   correct the signs
                # center3D    = x$signtile[[k]] * center3D    # recycling rule used here
                .Call( C_timesEqual, center3D, as.double( x$signtile[[k]] ), 2L )     # multiply in place

                idxsub  = idx[ masklarge ]

                #   assign only those centers for which deficit is LARGE
                centermax[ idxsub, ] = center3D[ masklarge,  ]
                }
            }


        out$pgramdf$centermax   = centermax

        out$pgramdf$shift   = sqrt( .rowSums( (out$pgramdf$centermax - out$pgramdf$center)^2 , length(betamax), 3 ) )

        deficient   =  tol <= out$pgramdf$shift     #  &  out$pgramdf$deficit != 0

        out$pgramdf$deficient   = deficient

        if( FALSE )
            {
            #   print some tracing data
            for( k in seq_len( length(x$zonogon) ) )
                {
                cat( "------------------- facet ", k, "  -----------------\n" )
                zono    = x$zonogon[[k]]

                subgnd  = getground(zono$matroid)

                idx = translateallpairs( subgnd, ground )

                print( out$pgramdf[idx, ] )
                }
            }


        #   compute area of all the pgrams
        crossprodsraw   = allcrossproducts( getmatrix(matsimple) )

        out$pgramdf$area    = sqrt( .colSums( crossprodsraw^2, nrow(crossprodsraw), ncol(crossprodsraw) ) )

        out$areadeficient   = 2*sum( out$pgramdf$area[deficient] )

        out$volumedeficient = (2/3) * sum( out$pgramdf$area[deficient] * out$pgramdf$deficit[deficient] )

        out$volume          = (2/3) * sum( out$pgramdf$area * betavec )

        if( FALSE )
            {
            mask    =  tol <= out$pgramdf$shift
            cat( "range of {shift < tol} =", range( out$pgramdf$shift[ ! mask ] ), "      (tol=", tol, ')\n' )
            cat( "range of {tol < shift} =", range( out$pgramdf$shift[ mask ] ), '\n' )
            }
        }

    return( out )
    }






#   inside2trans()
#
#   x           a zonohedron object
#   p           Mx3 matrix, with the M query points in the rowSums

#   value   a dataframe with columns
#       p               the given Mx3 input matrix
#       inside          TRUE means the linking number with the 2-transition surface is non-zero
#       linkingnumber   the linking number of the point w.r.t. the surface
#       distance        distance from the point to the surface
#       timecalc        time to do the calculation, in sec

#                       negative or 0 means inside, and positive means in the exterior
#                       NOTE: if positive then the distance is only approximate.

inside2trans <- function( x, p )       #  tol=5.e-12
    {
    if( ! inherits(x,"zonohedron") )
        {
        log_level( ERROR, "Argument x is invalid.  It is not a zonohedron." )
        return(NULL)
        }

    p   = prepareNxM( p, 3 )
    if( is.null(p) )    return(NULL)

    pgramdf   = allpgramcenters2trans( x )
    if( is.null(pgramdf) )  return(NULL)

    #   subtract x$center from all the given points
    #point_centered  = duplicate( p )
    #res = .Call( C_plusEqual, point_centered, -x$center, 1L )   # changes point_centered in place
    #if( is.null(res) )  return(NULL)
    point_centered  = .Call( C_sumMatVec, p, -x$center, 1L )

    matsimp = getsimplified( x$matroid )
    matgen  = getmatrix(matsimp)

    pointed = is_pointed(x)

    m   = nrow(p)

    inside      = rep( NA, m )      # logical
    linknum     = rep( NA_integer_, m )
    distance    = rep( NA_real_, m )
    timecalc    = rep( NA_real_, m )

    for( k in 1:m )
        {
        time_start  = gettime()
        
        pcent   = point_centered[k, ]

        #   quick check for black point or white point
        black_or_white = pointed  &&  ( all(pcent == -x$center)  ||  all(pcent == x$center) )

        if( ! black_or_white )
            #   the usual case
            distance[k] = .Call( C_dist2surface, matgen, pgramdf$idxpair, pgramdf$center, pgramdf$cross, pcent )
        else
            #   black and white points are vertices, special cases
            distance[k] = 0

        if( distance[k] != 0 )
            linknum[k] =   linkingnumber3( x, pgramdf, pcent )
            #   else if distance[k]==0 we just leave linknum[k] as it is, which is NA_integer_

        inside[k]   = (linknum[k] != 0L)

        timecalc[k] = gettime() - time_start
        }

    rnames  = rownames(p)
    if( is.null(rnames)  ||  anyDuplicated(rnames) )   rnames = 1:m

    out = data.frame( row.names=rnames )

    out$p               = p
    out$distance        = distance
    out$linkingnumber   = linknum
    out$inside          = inside
    out$timecalc        = timecalc

    return( out )
    }


#   x       a zonohedron object
#   type    'e' for edges, 'f' for facets, 'p' for points (centers of facets)
#   ecol    edge color
#   econc   if TRUE then concave edges overdrawn thick and red
#   fcol    color used for the coincident facets
#   falpha  opacity used for the coincident facets
#   normals if TRUE then facet normals are drawn
#   both    if TRUE draw both halves
#   bgcol   background color
#   add     if TRUE then add to an existing plot


plot2trans <- function( x, type='ef', ecol='black', econc=FALSE,
                                fcol='yellow', falpha=0.5, level=NULL,
                                normals=FALSE, both=TRUE, bgcol="gray40", add=FALSE, ... )
    {
    if( ! inherits(x,"zonohedron") )
        {
        log_level( ERROR, "Argument x is invalid.  It's not a zonohedron." )
        return(NULL)
        }

    if( ! requireNamespace( 'rgl', quietly=TRUE ) )
        {
        log_level( ERROR, "Package 'rgl' cannot be loaded. It is required for plotting the zonohedron." )
        return(FALSE)
        }


    matsimp = getsimplified( x$matroid )

    matgen  = getmatrix(matsimp)
    numgen  = ncol(matgen)

    if( ! is.null(level) )
        {
        #   check validity of level
        ok  = all( level %in% (0:(numgen-2L)) )
        if( ! ok )
            {
            log_level( ERROR, "argument level is invalid.  All values must be integers in [0,%d].", numgen-2 )
            return(FALSE)
            }
        }

    if( add )
        {
        if( rgl::cur3d() == 0 )
            {
            log_level( ERROR, "Cannot add surface to plot, because there is no rgl window open." )
            return(FALSE)
            }
        }
    else
        {
        #   start 3D drawing
        rgl::bg3d( color=bgcol )

        white   = 2 * x$center

        #cube    = rgl::scale3d( rgl::cube3d(col="white"), center[1], center[2], center[3] )
        #cube    = rgl::translate3d( cube, center[1], center[2], center[3]  )

        rgl::points3d( 0, 0, 0, col='black', size=10, point_antialias=TRUE )
        rgl::points3d( white[1], white[2], white[3], col='white', size=10, point_antialias=TRUE )
        rgl::points3d( x$center[1], x$center[2], x$center[3], col='gray50', size=10, point_antialias=TRUE )

        #   exact diagonal
        rgl::lines3d( c(0,white[1]), c(0,white[2]), c(0,white[3]), col=c('black','white'), lwd=3, lit=FALSE )
        }




    gndgen  = getground(matsimp)

    pgramdf = allpgramcenters2trans( x )

    #edgesok  = TRUE  # ! is.null(metrics$anglesDH)

    #if( grepl( 'e', type ) &&  ! edgesok )
    #    log.string( WARN, "Cannot draw edges because the dihedral angles are not available." )

    doedges = grepl( 'e', type )

    if( doedges && econc )      #&&  edgesok )
        {
        metrics = getmetrics2trans( x, tol=1.e-12 )
        if( is.null(metrics) )
            return(FALSE)

        anglesDH    = metrics$anglesDH

        #colvec  = ifelse( 0 <= anglesDH$angle, ecol, 'red' )
        #lwdvec  = ifelse( 0 <= anglesDH$angle, 1L, 3L )

        pivotmat = t( matgen[  , anglesDH$pivot ] )

        point0  = anglesDH$edgemid - 0.5*pivotmat
        point1  = anglesDH$edgemid + 0.5*pivotmat

        xyz = rbind( point0, point1 )

        m   = nrow(anglesDH)

        perm    = 1:(2*m)
        dim(perm)    = c(m,2L)
        perm = t(perm)
        dim(perm)    = NULL
        # print( perm )

        xyz = xyz[ perm, ]

        xyzdisp = .Call( C_sumMatVec, xyz, x$center, 1L )

        rgl::segments3d( xyzdisp, col=ecol, lwd=1 )

        conmask = (anglesDH$angle < 0)  #; print( conmask )

        if( econc && any(conmask) )
            {
            mask2   = rep(conmask,2)
            dim(mask2)  = c(m,2L)
            mask2   = t(mask2)
            dim(mask2)  = NULL

            rgl::segments3d( xyzdisp[mask2, ], col='red', lwd=3 )
            }

        if( both )
            {
            xyzdisp = .Call( C_sumMatVec, -xyz, x$center, 1L )

            #rgl::segments3d( xyzdisp, col=ecol, lwd=1 )

            if( econc && any(conmask) )
                rgl::segments3d( xyzdisp[mask2, ], col='red', lwd=3 )
            }
        }


    if( grepl( 'f', type ) )
        {
        #   draw filled quads

        pgrams  = nrow(pgramdf)
        step    = 4
        quadmat = matrix( 0, nrow=step*pgrams, ncol=3 )

        for( i in 1:pgrams )
            {
            center  = pgramdf$center[i, ]

            edge    = matgen[  , pgramdf$idxpair[i, ] ]  # 3x2 matrix

            k   = step*(i-1)

            quadmat[k+1, ] = center - 0.5 * edge[ , 1] - 0.5*edge[ , 2]
            quadmat[k+2, ] = center - 0.5 * edge[ , 1] + 0.5*edge[ , 2]
            quadmat[k+3, ] = center + 0.5 * edge[ , 1] + 0.5*edge[ , 2]
            quadmat[k+4, ] = center + 0.5 * edge[ , 1] - 0.5*edge[ , 2]
            }

        if( ! is.null(level) )
            {
            levelvec        = pgramdf$idxpair[ ,2] - pgramdf$idxpair[ ,1] - 1L
            #   repeat each value step (4) times
            levelvec        = matrix( levelvec, nrow=step, ncol=length(levelvec), byrow=TRUE )
            dim(levelvec)   = NULL      # back to vector
            levelmask       = levelvec %in% level
            levelmaskanti   = (numgen-2L - levelvec) %in% level
            }

        if( is.null(level) )
            xyz = .Call( C_sumMatVec, quadmat, x$center, 1L )
        else
            xyz = .Call( C_sumMatVec, quadmat[levelmask, ], x$center, 1L )

        #   compute the complementary color, and use this for the 'backfacing' facets
        fcol_comp   = grDevices::rgb( t( 255 - grDevices::col2rgb(fcol) ), max=255 )

        linknum = linkingnumber3( x, pgramdf )
        
        fcolvec = ifelse( 0 < linknum * pgramdf$beta, fcol, fcol_comp )

        #   replicate each color 4 times
        fcolmat = matrix( fcolvec, nrow=step, ncol=length(fcolvec), byrow=TRUE )
        fcolvec = as.character( fcolmat )
        
        if( ! is.null(level) )  fcolvec = fcolvec[levelmask]

        rgl::quads3d( xyz, col=fcolvec, alpha=falpha, lit=TRUE  )              # quad filled

        if( doedges )
            rgl::quads3d( xyz, col=ecol, lwd=1, front='lines', back='lines', lit=FALSE  )   # quad edges

        if( both )
            {
            if( is.null(level) )
                xyz = .Call( C_sumMatVec, -quadmat, x$center, 1L )
            else
                xyz = .Call( C_sumMatVec, -quadmat[levelmaskanti, ], x$center, 1L )

            rgl::quads3d( xyz, col=fcolvec, alpha=falpha, lit=TRUE )           # quad filled

            if( doedges )
                rgl::quads3d( xyz, col=ecol, lwd=1, front='lines', back='lines', lit=FALSE  )   # quad edges
            }
        }

    if( grepl( 'p', type ) )
        {
        #   draw centers
        xyz = .Call( C_sumMatVec, pgramdf$center, x$center, 1L )

        rgl::points3d( xyz[ ,1], xyz[ ,2], xyz[ ,3], col='black', size=3, point_antialias=TRUE )

        if( both )
            {
            xyz = .Call( C_sumMatVec, -pgramdf$center, x$center, 1L )
            rgl::points3d( xyz[ ,1], xyz[ ,2], xyz[ ,3], col='black', size=3, point_antialias=TRUE )
            }
        }

    if( normals )
        {
        xyz = .Call( C_sumMatVec, pgramdf$center, x$center, 1L )

        for( i in 1:nrow(xyz) )
            rgl::arrow3d( xyz[i, ], xyz[i, ] + pgramdf$cross[i, ], type="lines", col="black" )
        }

    return( invisible(TRUE) )
    }






#   arguments:
#
#   N       dimension of the cube, must be a positive integer
#   crange  range for the count of +1/2s for the edges, does not affect the vertex
#   type    'both' means both Type1 (BP)  and  Type2 (BS)
#           can also be 'BP'
#
#   returns list with components:
#       N       the input N
#       vertex  (N*(N-1)+2)x2 integer matrix with code for the vertex
#               the 1st int is the # of ones, and the 2nd is the starting position
#       edge    (2N*(N-2) + 2N) x 2 integer matrix with starting and ending index (in vertex) of the edges
#               this number applies only when crange=c(0L,N)
#

trans2subcomplex <- function( N, crange=c(0L,N), type='both' )
    {
    ok  = length(N)==1  &&  0<N
    if( ! ok )
        {
        log_level( ERROR, "N is invalid." )
        return(NULL)
        }

    ok  = is.numeric(crange)  &&  length(crange)==2  &&  0<=crange[1]  && crange[1]<crange[2]  && crange[2]<=N
    if( ! ok )
        {
        log_level( ERROR, "crange is invalid." )
        return(NULL)
        }

    vertex  = .Call( C_trans2vertex, N )
    colnames(vertex)    = c( "count", "start" )

    out = list()

    out$N       = N
    out$vertex  = vertex
    out$edge    = .Call( C_trans2edge, N, crange )

    if( type != 'both' )
        {
        rsums   = rowSums(vertex)

        if( toupper(type) == 'BP' )
            vvalid  = rsums <= N+1
        else if( toupper(type) == 'BS' )
            vvalid  = N+1 <= rsums  |  vertex[ ,2]==1       # add the vertices whose 1s start at position 1
        else
            {
            log_level( ERROR, "type=%s in invalid.", type )
            return(NULL)
            }

        vvalid = vvalid  |  is.na(rsums)   # is.na() is for the poles

        evalid  = vvalid[ out$edge[ ,1] ]  &  vvalid[ out$edge[ ,2] ]

        out$edge    = out$edge[ evalid, ]
        }

    #  print( out$edge )

    return(out)
    }



#   parameters:
#   N       the dimension of the cube, must be a positive integer
#   count   integer M-vector with the number of 1s in the vertex
#   start   integer M-vector with the starting index of the run of 1s, 1-based
#
#   returns an MxN matrix where the i'th row is the desired vertex of the N-cube, [-1/2,1/2]^N

vertexfromcode <- function( N, count, start )
    {
    M   = length(count)
    if( length(start) != M )    return(NULL)

    out = .Call( C_vertexfromcode, N, count, start )

    return(out)
    }

#   zono    a zonohedron, whose simplified matroid is generated by an 3 x N matrix of N generators,
#           defining a 2-transition surface in R^3
#
#   returns a data frame with data on the parallelograms in the 2-transition surface (not necessarily the boundary of zono)
#                   There are N*(N-1)/2 rows and these columns
#       idxpair     integer matrix with 2 columns i and j with 1 <= i < j <= n.  1-based
#       gndpair     idxpair converted to the ground set, another integer matrix with 2 columns
#       center      real matrix with 3 columns containing the corresponding pgram center,
#                   within the centered zonohedron.  These are computed efficiently using a cumulative matrix sum technique.
#       cross       a unit normal for the pgram, equal to the normalized cross-product of the 2 generators,
#                   and directly copied from the crossprods member of the matroid
#       beta        constant of the plane equation of the pgram.  These can be + or - or 0.
#
#       the row order is the same as the column order returned by allcrossproducts()

allpgramcenters2trans  <- function( zono )        #, centered=TRUE )
    {
    matsimple   = getsimplified(zono$matroid)

    matgen  = getmatrix(matsimple)

    #ok  = is.numeric(matgen)  &&  is.matrix(matgen)
    #if( ! ok )  return(NULL)

    n   = ncol(matgen)

    matcum  = .Call( C_cumsumMatrix, matgen, 2L )

    idxpair = .Call( C_allpairs, n )
    colnames(idxpair)   = c('i','j')

    #   pgrams  = nrow(idxpair)  #  N(N-1)/2

    #   center is loaded with the *uncentered* coords
    center  = .Call( C_allpgramcenters2trans, matgen, matcum )

    if( TRUE )
        {
        # translate from original to the *centered* zonohedron
        centerzono  = matcum[ ,n] / 2
        #center      = .Call( C_sumMatVec, center, -centerzono, 1L )

        .Call( C_plusEqual, center, -centerzono, 1L )       # change center[,] in place.  Might be a tiny bit faster

        #   print( str(center) )
        }

    out = data.frame( row.names=1:nrow(idxpair) )
    out$idxpair = idxpair

    ground  = getground( matsimple )
    gndpair         = ground[idxpair]
    dim(gndpair)    = dim(idxpair)
    out$gndpair     = gndpair

    out$center      = center

    # in the next line,  matsimple$crossprods is  3 x N(N-1)/2
    out$cross       = t(matsimple$crossprods)

    out$beta        = .rowSums( center * out$cross, nrow(center), ncol(center) )

    return( out )
    }


#   zono        the zonohedron
#   pgramdf     pgram data frame, as returned from allpgramcenters2trans(zono)
#   point       point from which to take the linking number, in centered zono coordinates
#               the default is the center of symmetry
#
#   returns the integer linking number.
#       if point[] is a vertex of the surface, it returns NA_integer_

linkingnumber <- function( zono, pgramdf=NULL, point=c(0,0,0) )
    {
    matsimp = getsimplified( zono$matroid )

    matgen  = getmatrix(matsimp)

    if( is.null(pgramdf) )
        {
        pgramdf   = allpgramcenters2trans( zono )
        if( is.null(pgramdf) )  return(NULL)
        }

    out = .Call( C_linkingnumber, matgen, pgramdf$idxpair, pgramdf$center, point )

    return( out )
    }

#   this version replaces pgramdf by the simpler matcum, but needs more work

linkingnumber2 <- function( zono, point=c(0,0,0) )
    {
    matsimp = getsimplified( zono$matroid )

    matgen  = getmatrix(matsimp)

    matcum  = cbind( 0, .Call( C_cumsumMatrix, matgen, 2L ) )

    out = .Call( C_linkingnumber2, matcum, point )

    return( out )
    }

#   this version replaces spherical triangles by planar polygons

linkingnumber3 <- function( zono, pgramdf=NULL, point=c(0,0,0) )
    {
    matsimp = getsimplified( zono$matroid )

    matgen  = getmatrix(matsimp)

    if( is.null(pgramdf) )
        {
        pgramdf   = allpgramcenters2trans( zono )
        if( is.null(pgramdf) )  return(NULL)
        }

    out = .Call( C_linkingnumber3, matgen, pgramdf$idxpair, pgramdf$center, point )

    return( out )
    }


#   k1 and k2   distinct integers in 1:n, not in any specific order
#   crossprods  3 x N*(N-1)/2 matrix of precomputed normalized cross products

crossprodlookup <- function( k1, k2, n, crossprods )
    {
    s   =  sign( k2 - k1 )
    if( s < 0 )
        {
        #   swap so that k1 < k2
        temp=k1 ; k1=k2 ; k2=temp
        }

    cp  = s * crossprods[ , (k1-1)*n - ((k1)*(k1+1))/2 + k2  ]    #; cat( "cp=", cp )

    return( cp )
    }




#   k0          index of the pivot edge
#   k1, sign1   index and sign of wing #1
#   k2, sign2   index and sign of wing #2
#   matgen      3 x N matrix of edge generators
#   crossprods  3 x N*(N-1)/2 matrix of precomputed normalized cross products

hingeangle  <- function( k0, k1, sign1, k2, sign2, matgen, crossprods )
    {
    n   = ncol(matgen)

    signwings   = sign1 * sign2

    cp1 = signwings * crossprodlookup( k1, k0, n, crossprods )      #    s * crossprods[ , PAIRINDEX( k1, k0, n ) ]

    cp2 = crossprodlookup( k0, k2, n, crossprods )      #s * crossprods[ , PAIRINDEX( k0, k2, n ) ]

    theta   = angleBetween( cp1, cp2, unitized=TRUE )   #; cat( "theta=", theta, '\n' )

    cp  = signwings * crossprodlookup( k1, k2, n, crossprods )      # s * crossprods[ , PAIRINDEX( k1, k2, n ) ]    ; cat( "cp=", cp )

    s   = sign( sum( matgen[ ,k0] * cp ) )              #    ; cat( "  gen=", matgen[ ,k],  "    s=", s, '\n' )

    return( s * theta )
    }

#   zono        a zonohedron
#   signsurf    if k1<k2 then the crossprod lookup is outward when signsurf=1 and inward when signsurf=-1

#   returns data.frame with N*(N-1)  rows and these columns
#       pivot   integer index of the generator where dihedral angle pivots, the pivot of the "hinge"
#       wing1   index of the generator forming wing 1 of the "hinge"
#       wing2   index of the generator forming wing 2 of the "hinge"
#       level   the # of 1s in the level where the edge ends
#       angle   the external dihedral angle, where positive is convex and negative is concave
#       edgemid midpoint of the edge, in the *centered* polyhedron


edgeangles2trans  <- function( zono, signsurf )
    {
    matsimp     = getsimplified( zono$matroid )
    matgen      = getmatrix( matsimp )
    crossprods  = matsimp$crossprods

    ok  = is.numeric(matgen)  &&  is.matrix(matgen)  &&  nrow(matgen)==3
    if( ! ok )  return(NULL)

    gensum  = .Call( C_cumsumMatrix, matgen, 2L )

    n   = ncol(matgen)

    knext   = c( 2L:n, 1L )
    kprev   = c( n, 1L:(n-1L) )

    hinges  =  n*(n-1L)     # later n

    pivot   = rep( NA_integer_, hinges )
    wing1   = matrix( NA_integer_, hinges, 2 )  ; colnames(wing1) = c( "index", "sign" )
    wing2   = matrix( NA_integer_, hinges, 2 )  ; colnames(wing2) = colnames(wing1)
    level   = rep( NA_integer_, hinges )
    angle   = rep( NA_real_, hinges )
    edgemid = matrix( NA_real_, hinges, 3 )

    #   group # is for the bottom level 0, for black and white points
    for( k in 1:n )
        {
        pivot[k]    = k
        k1          = kprev[k]
        k2          = knext[k]

        wing1[k, ]  = c(k1,+1L)
        wing2[k, ]  = c(k2,+1L)

        level[k]    = 0L

        angle[k] = hingeangle( k, k1, +1, k2, +1, matgen, crossprods )

        edgemid[k, ]    = 0.5 * matgen[ ,k]
        }


    #   group #2 is for the higher levels, away from the black point
    kwrap   = rep( 1:n, 2 )
    white   = gensum[ ,n]   #;  cat( "white=", white, '\n' )

    count   = n
    for( shift in 1:(n-2) )
        {
        for( k in 1:n )
            {
            k1  = kwrap[ k+shift ]
            k2  = kwrap[ k1+1 ]

            count   = count+1L

            pivot[count]    = k
            wing1[count, ]  = c(k1,-1L)
            wing2[count, ]  = c(k2,+1L)

            #cat( "k=", k, "   k2=", k2, '\n' )
            mid = 0.5 * matgen[ ,k]

            level[count]    = shift

            if( k+shift <= n )
                {
                mid = mid  +  (gensum[ ,k1] - gensum[ ,k])
                }
            else
                {
                #   wrapped around
                mid = mid  +  (gensum[ ,k-1L] - gensum[ ,k1])
                mid = white - mid
                }

            angle[count]    = hingeangle( k, k1, -1, k2, +1, matgen, crossprods )

            #cat( "mid=", mid, '\n' )
            edgemid[count, ]    = mid
            }
        }

    out = data.frame( row.names=1:hinges )
    out$pivot   = pivot
    out$wing1   = wing1
    out$wing2   = wing2
    out$level   = level
    out$angle   = signsurf * angle
    out$edgemid = .Call( C_sumMatVec, edgemid, -white/2, 1L ) # translate to the *centered* polyhedron

    return( out )
    }

#   x       a zonohedron
#   hpidx   index of a non-trivial hyperplane in x
#
#   returns TRUE or FALSE

is_2transfacetNT  <- function( x, hpidx )
    {
    matsimple   = getsimplified(x$matroid)

    numgen  = length( matsimple$hyperplane[[hpidx]] )

    if( numgen <= 2 )
        {
        log_level( WARN, "internal error.  hpidx=%d is a trivial %d-point hyperplane.", hpidx, numgen )
        return( FALSE )
        }

    zgon    = x$zonogon[[hpidx]]
    if( ! is_salient(zgon) )    return(FALSE)


    #   check that the generators are monotone ordered by angle
    mat = getmatrix( getsimplified(zgon$matroid) )

    theta   = atan2( mat[2, ], mat[1, ] )

    #   rotate so facet0 is at theta=0
    theta   = theta - theta[ zgon$facet0[1] ]   #; print( theta )

    #   wrap to range[-pi,pi]
    theta   = ((theta+pi) %% (2*pi)) - pi       #; print( theta )

    perm        = order( theta )
    n           = length(perm)
    monotone    = all(perm==1:n) || all(perm==n:1)
    if( ! monotone )    return(FALSE)


    #   check that the generators of the zonogon are contiguous in those of the zonohedron, with wrap-around (cyclic)
    ground      = getground(matsimple)  # strictly increasing

    subground   = getground( getsimplified(zgon$matroid) )

    if( ! is_contiguous( subground, ground ) )  return( FALSE )

    return( TRUE )
    }


is_contiguous   <- function( subground, ground, cyclic=TRUE )
    {
    idx = match( subground, ground )    #; print(idx)

    if( any( is.na(idx) ) )
        {
        log_level( WARN, "internal error. ground set of non-trivial facet is not a subset of the zonohedron ground set." )
        return( FALSE )
        }

    diffidx = diff( idx )

    count1  = sum( diffidx == 1 )

    # m   = length(ground)
    n   = length(subground)

    if( cyclic )
        out = (count1 == n-1)  ||  ( (count1 == n-2) && (sum(diffidx==(1-length(ground))) == 1) )
    else
        out = (count1 == n-1)

    return( out )
    }




#   idxpair     pair of distinct integer indexes, between 1 and n.  NOT verified
#   alpha       pair of points in [0,1].  NOT verified
#   n           integer n >= 3.  NOT verified
#
#   returns 2-transition point in n-cube with the given data

pcubefromdata   <- function( idxpair, alpha, n )
    {
    out = numeric( n )

    out[idxpair]  = alpha

    if( idxpair[1] < idxpair[2] )
        {
        #   Type I
        if( idxpair[1]+1 <= idxpair[2]-1  ) out[ (idxpair[1]+1):(idxpair[2]-1) ]    = 1
        }
    else
        {
        #   Type II
        if( idxpair[1]+1 <= n ) out[ (idxpair[1]+1):n ] = 1

        if( 1 <= idxpair[2]-1 ) out[ 1:(idxpair[2]-1) ] = 1
        }

    return( out )
    }


#   x           a zonohedron object
#   base        basepoint of all the rays, a 3-vector
#   direction   M x 3 matrix with non-zero directions in the rows
#   invert      return 2-transition point in the cube
#   plot        add 3D plot
#   tol         tolerance for being strictly starshaped, and for intersection with a pole
#
#   value   a dataframe with columns
#           base        given basepoint of all the rays (all the same)
#           direction   given directions of the rays
#           gndpair     the 2 indexes of the generators of the pgram that the ray intersects
#           alpha       the 2 coordinates of the intersection point within the pgram
#           tmax        ray parameter of intersection with pgram
#           point       point of intersection of the ray and the pgram
#           iters       the number of pgrams searched, until the desired one was found
#           timetrace   time to do the raytrace, in seconds
#
#   and if invert is TRUE, then these:
#           pcube       a point in the cube that maps to the point on the surface
#

raytrace2trans <- function( x, base, direction, invert=FALSE, plot=FALSE, tol=1.e-12, ... )
    {
    if( ! inherits(x,"zonohedron") )
        {
        log_level( ERROR, "Argument x is invalid.  It's not a zonohedron." )
        return(NULL)
        }

    ok  = is.numeric(base)  &&  length(base)==3  &&  all( is.finite(base) )
    if( ! ok )
        {
        log_level( ERROR, "base is invalid. It must be a numeric vector of length 3, and all entries finite." )
        return(NULL)
        }

    center  = getcenter(x)

    base_centered = base - center

    base_at_pole    = all(base_centered == -center)  ||  all(base_centered == center)
    if( base_at_pole )
        {
        log_level( ERROR, "The base=(%g,%g,%g) of the rays is at a pole, which cannot be processed at this time.",
                            base[1], base[2], base[3] )
        return(NULL)
        }

    direction   = prepareNxM( direction, 3 )
    if( is.null(direction) )    return(NULL)

    if( any( x$matroid$multiplesupp$mixed ) )
        {
        log_level( ERROR, "In the zonohedron generators, one of the multiple groups has mixed directions.  The 2-transition surface cannot be processed in this version." )
        return(NULL)
        }


    symmetric   = all( base_centered == 0 )

    pgramdf = allpgramcenters2trans(x)
    if( is.null(pgramdf) )  return(NULL)

    linknum = linkingnumber3( x, pgramdf, base_centered )

    if( is.na(linknum) )
        {
        log_level( ERROR, "The linking number is undefined for point=(%g,%g,%g).",
                            base[1], base[2], base[3] )
        return(NULL)
        }

    if( abs(linknum) != 1 )
        {
        log_level( ERROR, "The linking number at basepoint=(%g,%g,%g) is %d, which is not +1 or -1.  The ray intersection is never unique.",
                                base[1], base[2], base[3], linknum )
        return(NULL)
        }




    matsimple   = getsimplified( x$matroid )


    if( ! all( x$matroid$multiplesupp$contiguous ) )
        {
        log_level( WARN, "The 2-transition surface is not starshaped at any point, because one of the multiple groups is not contiguous. The ray intersection may not be unique." )
        }
    else
        {
        if( linknum == -1 )
            {
            #   use antipodal facets instead
            #   so we can compare beta with the corresponding betamax from the zonohedron,
            #   we only do this when starshaped is TRUE, see below
            #cat( "linknum =", linknum, "  so reversing beta.\n" )
            #   pgramdf$center = -(pgramdf$center)
            pgramdf$beta   = -(pgramdf$beta)
            }

        betamin = min( pgramdf$beta )
        if( betamin < tol )
            {
            log_level( WARN, "The 2-transition surface is not strictly starshaped at any point, because min(beta) = %g < %g. The ray intersection may not be unique.",
                                betamin, tol )
            #   return(NULL)
            }
        else
            {
            basedotnormal   = as.double( base_centered  %*%  matsimple$crossprods )  #; print( str(basedotnormal) )

            ok  = all( abs(basedotnormal) < pgramdf$beta - tol )
            if( ! ok )
                {
                log_level( WARN, "The 2-transition surface is not strictly starshaped at point=(%g,%g,%g). The ray intersection may not be unique.",
                                     base[1], base[2], base[3] )
                #return(NULL)
                }
            }
        }


    matgen  = getmatrix( matsimple )

    m       = ncol(matgen)

    # extend pgramdf with antipodal data
    #pgramdf     = pgramdf_plus( pgramdf, base_centered )

    #print( str(pgramdf) )

    #   extend pairs in usual i<j order, with new pairs where i>j
    #   the number of rows is m*(m-1)
    idxpair_plus    = rbind( pgramdf$idxpair, pgramdf$idxpair[ ,2:1] )

    rays    = nrow(direction)

    matcum  = cbind( 0, .Call( C_cumsumMatrix, matgen, 2L ) )       # this has m+1 columns

    # print( gensum )

    tmax        = rep(NA_real_,rays)
    idxpair     = matrix(NA_integer_,rays,2)
    alpha       = matrix(NA_real_,rays,2)
    point       = matrix(NA_real_,rays,3)
    iters       = rep( NA_integer_, rays )
    transitions = rep( NA_integer_, rays )
    timetrace   = rep(NA_real_,rays)

    if( invert )    pcube = matrix( NA_real_, rays, m )

    for( k in 1:rays )
        {
        #   cat( "k =", k, '\n' ) ; flush.console()

        time_start  = gettime()

        #print( matgen[ ,cw_init[1], drop=F ] )

        #   get current direction
        dir = direction[k, ]

        if( all(dir==0) )   next    # ray is undefined

        dat = raytrace2trans_single( base, dir, center, matgen, matcum, pgramdf, idxpair_plus, tol=tol )

        timetrace[k]    = gettime() - time_start

        # print( str(dat) )

        if( ! is.null(dat) )
            {
            idxpair[k, ]    = dat$idxpair
            alpha[k, ]      = dat$alpha
            tmax[k]         = dat$tmax
            point[k, ]      = dat$point
            iters[k]        = dat$iters

            transitions[k]  = 2L

            if( invert )    pcube[k, ]  = pcubefromdata( dat$idxpair, dat$alpha, m )
            }

        #print( (dat$XYZ - base) / rtdata$direction[k, ] )
        }


    rnames  = rownames(direction)
    if( is.null(rnames)  ||  anyDuplicated(rnames) )   rnames = 1:rays

    #   convert idxpair to gndpair, for output
    gndpair         = getground( matsimple )[ idxpair ]
    dim(gndpair)    = dim(idxpair)

    out = data.frame( row.names=rnames )

    out$base        = matrix( base, rays, 3, byrow=TRUE )  # replicate base to all rows
    out$direction   = direction
    out$gndpair     = gndpair
    #out$idxpair     = idxpair
    out$alpha       = alpha
    out$tmax        = tmax
    out$point       = point
    #out$transitions = transitions
    out$iters       = iters
    out$timetrace   = timetrace

    if( invert )    out$pcube   = invertcubepoints( x, pcube )

    if( plot )
        {
        if( ! requireNamespace( 'rgl', quietly=TRUE ) )
            log_level( WARN, "Package 'rgl' is required for plotting.  Please install it." )
        else if( rgl::cur3d() == 0 )
            log_level( WARN, "Cannot add raytrace to plot, because there is no rgl window open." )
        else
            {
            xyz = matrix( base, nrow=rays, ncol=3, byrow=TRUE )
            xyz = rbind( xyz, point )

            perm    = 1:(2*rays)
            dim(perm)    = c(rays,2L)
            perm = t(perm)
            dim(perm)    = NULL
            # print( perm )

            xyz = xyz[ perm, ]

            col = 'red'

            rgl::segments3d( xyz[ ,1], xyz[ ,2], xyz[ ,3], col=col )

            rgl::points3d( base[1], base[2], base[3], col=col, size=5, point_antialias=TRUE )
            rgl::points3d( point[ ,1], point[ ,2], point[ ,3], col=col, size=5, point_antialias=TRUE )

            out = invisible(out)
            }
        }

    return(out)
    }

#   base            base of ray, uncentered
#   dir             direction of ray
#   center          center of of symmetry
#   matgen          3 x M matrix of generators
#   matcum          3 x M+1 matrix = cumsum of matgen
#   pgramdf         data frame with idxpair, center, and other variables.   not extended with antipodal data
#   idxpair_plus    integer matrix with 2 columns, extended with antipodal data
#   tol             tolerance for poletest()

#   returns a list with these items:
#       idxpair     indexes of the generators of the pgram, or both NA if the ray intersects a pole
#       alpha       2 coords in [0,1] for point inside the pgram
#       iters       # of iterations, or 0 when the ray intersect a pgram that intersects a pole
#       tmax        parameter of ray when it intersects the surface, or NA in case of failure
#       point       intersection of ray with the surface, or NA in case of failure

raytrace2trans_single   <- function( base, dir, center, matgen, matcum, pgramdf, idxpair_plus, tol=1.e-12 )
    {
    m   = ncol(matgen)

    base_centered   = base - center

    #   from current dir, compute a 2x3 projection matrix
    frame3x3    = goodframe3x3( dir )

    matproj     = t( frame3x3[ , 1:2 ] )

    success = FALSE

    iters   = 0L

    #cat( "------------   dir =", dir, "     -------------", '\n' )


    ########   Phase 1 - test for intersection at or near a pole

    # dat = poletest( base_centered, dir, center, matgen, matproj, tol=tol )

    dat = poletest_v2( base_centered, dir, center, matproj, tol=tol )

    if( ! is.null(dat) )
        {
        #   test for a good intersection at or near a pole
        alpha   = dat$alpha

        if( is.na( dat$idxpair[1] ) )
            {
            #   intersection at a pole
            dat$XYZ = (1 +  dat$signpole) * center   # so either 0 or 2*center
            success = TRUE
            }
        else
            {
            #   test for intersection with a pgram that intersects a pole
            ok  = all( 0 <= alpha )  &&  all( alpha <= 1 )

            if( ok )
                {
                #   compute XYZ taking advantage of all the 0s
                mat3x2  = matgen[ , dat$idxpair ]

                if( dat$signpole == -1 )
                    #   near 0
                    dat$XYZ = as.double( mat3x2 %*% alpha )
                else
                    #   near white
                    dat$XYZ = 2*center - as.double( mat3x2 %*% (1-alpha) )    # complement

                #   final check on orientation
                #   (dat$XYZ - base) must be parallel to dir
                if( 0 < sum( (dat$XYZ - base)*dir ) )
                    success = TRUE
                }
            }

        if( success )
            {
            out = list()

            out$idxpair = dat$idxpair
            out$alpha   = alpha
            out$iters   = 0L

            # redo to match documentation
            out$tmax    = sum( (dat$XYZ - base)*dir ) / sum( dir^2 )
            out$point   = base  +  out$tmax * dir

            return( out )
            }
        }


    if( FALSE )
        {
        #   use BRUTE FORCE search in 3D
        #   only useful for timing comparison
        res = findpgram3D( base, dir, center, matgen, pgramdf )

        if( ! is.null(res) )
            {
            # print( str(res) )   ; flush.console()

            out = list()

            out$idxpair = res$idxpair
            out$alpha   = res$alpha
            out$iters   = iters + res$idx

            #   an alternate way to get XYZ
            #   no antipodal data, so res$idx <= nrow(pgramdf)
            centerpg    = pgramdf$center[ res$idx, ]
            XYZ     = center + centerpg + as.double( matgen[ , res$idxpair ] %*% (res$alpha - 1/2) )

            #   typical way to get XYZ
            # XYZ = XYZfrom2trans( res$idxpair, res$alpha, matgen, matcum )  # ; print( XYZ )

            # redo to match documentation
            out$tmax    = sum( (XYZ - base)*dir ) / sum( dir^2 ) # sqrt( sum( (XYZ - base)^2 ) / sum( dir^2 ) )
            out$point   = base  +  out$tmax * dir

            return( out )
            }
        }


    #   a near-pole intersection did not work
    #   rotate face centers to align direction with z-axis

    temp    = pgramdf$center %*% frame3x3

    #   extend with antipodal centers, the number of rows is now m*(m-1)
    #time_bind   = gettime()

    #centerrot   = rbind( temp, -temp )

    centerrot   = .Call( C_extend_antipodal, temp )

    #cat( "bind time =", gettime() - time_bind, '\n' )


    baserot     = as.double( base_centered %*% frame3x3 )
    # cat( "baserot = ", baserot, '\n' )

    if( FALSE )
        {
        #   find a good initial pgram2 for iteration
        time_initial   = gettime()

        if( FALSE )
            {
            qpoint  = baserot[1:2]

            matdiff     = pgram2df$center - matrix( qpoint, nrow(pgram2df), 2, byrow=TRUE )

            test = .rowSums( matdiff^2, nrow(matdiff), 2L )

            #   but force skipping of unit vectors in the wrong halfspace
            test[ pgram2df$z < 0 ]    = Inf

            imin    = which.min(test)

            #cat( "1st imin =", imin, '\n' )
            }

        #   ignore facets whose centers have negative z coords
        baserot[3]  = max( baserot[3], 0 )

        #   call optimized C function to get the closest center to qpoint
        imin = .Call( C_optimalcenter, centerrot, baserot ) + 1L   # 0-based to 1-based

        if( length(imin)==0 )
            {
            #   C_optimalcenter() failed !
            #   all the centers are below baserot[3], which is unusual
            #   there must be a few large faces, e.g. a cube
            #   choose the center with largest z
            imin    = which.max( centerrot[ ,3] )

            #cat( "center with largest z:  imin =", imin, "    initial idxpair0 =", idxpair_plus[imin, ], '\n' )
            }

        #cat( "imin =", imin, '\n' )

        idxpair0    = idxpair_plus[imin, ]

        #cat( "idxpair0 =", idxpair0, "  centermin =", centerrot[imin, ], '\n' )

        #cat( "Initial time =", gettime() - time_initial, '\n' )

        #####    use custom iteration on 2D pgrams, to find the one that contains the query point

        matgen2 = matproj %*% matgen

        res = findpgram2D( centerrot, baserot, matgen2, idxpair0 )   # idxpair0

        if( ! is.null(res) )
            {
            out = list()

            out$idxpair = res$idxpair
            out$alpha   = res$alpha
            out$iters   = iters + res$iters

            XYZ = XYZfrom2trans( res$idxpair, res$alpha, matgen, matcum )  # ; print( XYZ )

            #XYZ = XYZfrom2trans_slow( res$idxpair, res$alpha, matgen ) ; print( XYZ )

            # redo to match documentation
            out$tmax    = sum( (XYZ - base)*dir ) / sum( dir^2 ) # sqrt( sum( (XYZ - base)^2 ) / sum( dir^2 ) )
            out$point   = base  +  out$tmax * dir

            return( out )
            }
        }

    #  res = findpgram3D( base, dir, center, matgen, pgramdf )


    if( TRUE )
        {
        #####   use brute force search over all suitable pgrams,
        #       but using C function for great speed

        #   ignore facets whose centers have negative z coords
        #   baserot[3]  = max( baserot[3], 0 )

        genrot = crossprod( frame3x3, matgen )      # same as t(frame3x3) %*% matgen

        res = findpgram2D_v3( centerrot, baserot, idxpair_plus, genrot )

        if( ! is.null(res) )
            {
            out = list()

            out$idxpair = res$idxpair
            out$alpha   = res$alpha
            out$iters   = iters + res$idx

            XYZ = XYZfrom2trans( res$idxpair, res$alpha, matgen, matcum )  # ; print( XYZ )

            #XYZ = XYZfrom2trans_slow( res$idxpair, res$alpha, matgen ) ; print( XYZ )

            # redo to match documentation
            out$tmax    = sum( (XYZ - base)*dir ) / sum( dir^2 ) # sqrt( sum( (XYZ - base)^2 ) / sum( dir^2 ) )
            out$point   = base  +  out$tmax * dir

            return( out )
            }

        }


    #  res = findpgram3D( base, dir, center, matgen, pgramdf )

    return( NULL )
    }





XYZfrom2trans   <- function( idxpair, alpha, matgen, matcum )
    {
    m   = ncol(matgen)

    slo     = (idxpair[1]    - alpha[1]) %% m
    shi     = (idxpair[2]-1  + alpha[2]) %% m

    ilo = as.integer( slo + 1 )
    ihi = as.integer( shi + 1 )

    XYZlo   = matcum[ , ilo ] + (slo-floor(slo))*matgen[ , ilo ]
    XYZhi   = matcum[ , ihi ] + (shi-floor(shi))*matgen[ , ihi ]

    XYZ = as.double( XYZhi - XYZlo )

    if( ihi < ilo )
        # this is a bandstop, Type II
        XYZ = XYZ + matcum[ ,m+1]

    return( XYZ )
    }

XYZfrom2trans_slow  <- function( idxpair, alpha, matgen )
    {
    pcube   = pcubefromdata( idxpair, alpha, ncol(matgen) )

    out     = as.double( matgen %*% pcube )

    return( out )
    }



#   base    basepoint of the ray, centered
#   dir     direction of the ray
#   pole    the "positive" pole, with all coeffs +1/2, the "negative" pole is -pole, centered
#   matgen  3 x M matrix of generators
#   matproj 2 x 3 projection matrix to plane normal to dir
#   tol     tolerance for the ray going through a pole
#
#   first tests for ray passing within tol of a pole
#   next finds which one of the M sectors at that pole that the ray intersects
#       the 2 alpha coefficients of that intersection point are non-negative

poletest <- function( base, dir, pole, matgen, matproj, tol=0 )
    {
    time_start  = gettime()

    #cat( "base =", base, "   dir =", dir, "   +pole =", pole, '\n' )

    pd  = sum( pole*dir )
    bd  = sum( base*dir )

    if( abs(pd) <= bd )
        {
        #   both pole and -pole are on the negative side of the projection plane
        #cat( "both poles on - side of plane.", '\n' )
        return(NULL)
        }


    #   determine signpole -- select either +pole or -pole
    if( pd <= bd )
        #   pole is on the negative side, so -pole is on the positive side
        signpole    = -1L
    else if( -pd <= bd )
        #   -pole is on the negative side, so +pole is on the positive side
        signpole    = +1L
    else
        {
        #   both -pole and +pole are on the positive side
        #   choose the closest pole to base, after projection
        poleproj    = as.double( matproj %*% pole )
        baseproj    = as.double( matproj %*% base )

        signpole = sign( sum(poleproj*baseproj) )
        if( signpole == 0 ) signpole = -1

        #cat( "both poles on + side of plane.  signpos =", signpole, '\n' )
        }

    #dirpole = sum( dir * pole )

    #dirbase = sum( dir * base )
    #if( dirpole < dirbase ) return(NULL)    # pole is not on + side of the hyperplane

    #   get the query point
    #   in this case we choose a translation so the pole is at 0
    q   = as.double( matproj %*% (base - signpole*pole) )

    normq   = sum( abs(q) )   # ;

    #cat( "signpole =", signpole,  "  q =", q, "normq =", normq, "    tol =", tol, '\n' )

    if( normq <= tol )
        {
        #   ray passes right through one of the poles
        out = list()
        out$signpole    = signpole
        out$idxpair     = c(NA_integer_,NA_integer_)
        out$alpha       = rep( (signpole+1)/2, 2 )      # -1 -> 0    and   +1 -> 1
        return( out )
        }


    #   project generators onto plane, and translate
    #   at the positive pole (signpole==1) we want to *reverse* direction
    gen2    = (-signpole * matproj) %*% matgen     # gen2 is 2 x M

    #   find the norm of all these 2D generators
    normgen = .colSums( abs(gen2), nrow(gen2), ncol(gen2) )

    if( any( normgen == 0 ) )
        {
        log_level( FATAL, "Internal error.  %d of the %d generators project to 0.",
                            sum(normgen==0), length(normgen) )
        return( NULL )
        }

    #cat( "normq = ", normq, "   max(normgen) =", max(normgen), '\n' )

    if( 2 * max(normgen) < normq )
        #   q is too large
        return(NULL)

    thetavec    = atan2( gen2[2, ], gen2[1, ] )
    #cat( "thetavec =", thetavec, '\n' )


    thetadiff   = diff(thetavec)
    #cat( "thetadiff.   neg:", sum(thetadiff<0), "  pos:", sum(0<thetadiff), "  zero:", sum(thetadiff==0), '\n' )


    theta   = atan2( q[2], q[1] )

    #   find the generator closest to q in angle
    imin    = which.min( abs(theta - thetavec) )

    #cat( "theta = ", theta, "   imin =", imin,  "  thetavec[imin] =", thetavec[imin], '\n' )

    m   = length(thetavec)

    iprev = ((imin-2L) %% m ) + 1L
    inext = ( imin %% m ) + 1L

    if( on_arc( theta, c(thetavec[iprev],thetavec[imin]) ) )
        iother = iprev
    else if( on_arc( theta, c(thetavec[inext],thetavec[imin]) ) )
        iother = inext
    else
        {
        log_level( FATAL, "Internal error. Cannot find angular interval containing theta=%g.", theta )
        return(NULL)
        }

    idxpair = c(imin,iother)


    #   get the order right
    if( 0 < signpole )
        {
        #   ray nearly intersects positive pole, so cube point is mostly 1s
        #   idxpair should be decreasing, except when 1,m
        ok  = idxpair[2] < idxpair[1]

        if( all( idxpair %in% c(1,m) ) )    ok = ! ok

        if( ! ok )  idxpair = idxpair[ 2:1 ]    # swap
        }
    else
        {
        #   ray nearly intersects negative pole, so cube point is mostly 0s
        #   idxpair should be increasing, except when m,1
        ok  = idxpair[1] < idxpair[2]

        if( all( idxpair %in% c(m,1) ) )    ok = ! ok

        if( ! ok )  idxpair = idxpair[ 2:1 ]    # swap
        }


    log_level( TRACE, "Found interval idxpair=%d,%d.  signpole=%g\n", idxpair[1], idxpair[2], signpole )


    mat2x2  = gen2[ , idxpair ]  #; print( mat2x2 )

    alpha   = as.double( solve( mat2x2, q ) )

    ok  = all( 0 <= alpha )

    #cat( "alpha =", alpha, "  ok =", ok, '\n' )

    if( ! ok )
        {
        log_level( FATAL, "Internal error.  alpha=%g,%g and one is negative.", alpha[1], alpha[2] )
        return( NULL )
        }

    out = list()
    out$signpole    = signpole
    out$idxpair     = idxpair
    out$alpha       = alpha

    if( signpole == 1 ) out$alpha = 1 - out$alpha

    #cat( "poletest(). time =", gettime() - time_start, '\n' )

    return( out )
    }



#   base    basepoint of the ray, centered
#   dir     direction of the ray
#   pole    the "positive" pole, with all coeffs +1/2, the "negative" pole is -pole, centered
#   matproj 2 x 3 projection matrix to plane normal to dir
#   tol     tolerance for the ray going through a pole
#
#   This is v. 2.   It only tests for ray passing within tol of a pole

poletest_v2 <- function( base, dir, pole, matproj, tol=0 )
    {
    time_start  = gettime()

    #cat( "base =", base, "   dir =", dir, "   +pole =", pole, '\n' )

    pd  = sum( pole*dir )
    bd  = sum( base*dir )

    if( abs(pd) <= bd )
        {
        #   both pole and -pole are on the negative side of the projection plane
        #cat( "both poles on - side of plane.", '\n' )
        return(NULL)
        }


    #   determine signpole -- select either +pole or -pole
    if( pd <= bd )
        #   pole is on the negative side, so -pole is on the positive side
        signpole    = -1L
    else if( -pd <= bd )
        #   -pole is on the negative side, so +pole is on the positive side
        signpole    = +1L
    else
        {
        #   both -pole and +pole are on the positive side
        #   choose the closest pole to base, after projection
        poleproj    = as.double( matproj %*% pole )
        baseproj    = as.double( matproj %*% base )

        signpole = sign( sum(poleproj*baseproj) )
        if( signpole == 0 ) signpole = -1

        #cat( "both poles on + side of plane.  signpos =", signpole, '\n' )
        }

    #   get the query point
    #   in this case we choose a translation so the pole is at 0
    q   = as.double( matproj %*% (base - signpole*pole) )

    normq   = sum( abs(q) )   # ;

    #cat( "signpole =", signpole,  "  q =", q, "normq =", normq, "    tol =", tol, '\n' )

    if( normq <= tol )
        {
        #   ray passes right through one of the poles
        out = list()
        out$signpole    = signpole
        out$idxpair     = c(NA_integer_,NA_integer_)
        out$alpha       = rep( (signpole+1)/2, 2 )      # -1 -> 0    and   +1 -> 1
        return( out )
        }

    #   failed to intersect the pole
    return( NULL )
    }


pgramdf_plus    <- function( pgramdf, base_centered=c(0,0,0) )
    {
    # extend pgramdf with antipodal data
    out = data.frame( row.names=1:nrow(pgramdf) )
    out$idxpair   = pgramdf$idxpair[ , 2:1]     # swap
    out$gndpair   = pgramdf$gndpair[ , 2:1]     # swap
    out$center    = -(pgramdf$center)
    out$beta      = -(pgramdf$beta)

    out = rbind( pgramdf, out )

    if( FALSE )
        {
        #   subtract base_centered from all pgram centers
        #   these points on the unit sphere will be used later
        #xyz = duplicate( out$center )
        #res = .Call( C_plusEqual, xyz, -base_centered, 1L )   # changes xyz in place
        #if( is.null(res) )  return(NULL)

        xyz = .Call( C_sumMatVec, out$center, -base_centered, 1L )

        #   and then unitize
        ok  = .Call( C_normalizeMatrix, xyz, 1L )   # changes xyz in place
        if( ! ok )  return(NULL)

        #   add xyz as a new column, to be used later
        out$unit    = xyz
        }


    return( out )
    }

#   given a tiling by pgrams of a part of the plane, and a query point
#   find the pgram that contains the query point
#
#   centerrot   N*(N-1) x 3 matrix with pgram centers,
#               only the first 2 columns are used, the 3rd z-coordinate might be used in future testing
#   baserot     3-vector with query point, 3rd z-coord not used
#               the goal is to find the pgram that contains this point
#   matgen2     2 x N matrix of generators
#   idxpair0    integer pair to start the search

findpgram2D   <- function( centerrot, baserot, matgen2, idxpair0 )
    {
    time_start = gettime()

    n   = ncol(matgen2)

    if( nrow(centerrot) != n*(n-1) )
        {
        log_level( ERROR, "nrow(centerrot) = %g != %g = n*(n-1).  n=%g", nrow(centerrot), n*(n-1), n )
        return(NULL)
        }

    qpoint  = baserot[1:2]  #;       cat( "qpoint = ", qpoint, '\n' )

    #   use variable i and j to abbreviate 2 ints in idxpair0
    i   = idxpair0[1]
    j   = idxpair0[2]

    maxiters    = as.integer( max( 0.75*n, 50 ) )
    success     = FALSE

    #   make increment and decrement vectors
    idx_inc = c(2:n,1)
    idx_dec = c(n,1:(n-1))

    for( iter in 1:maxiters )
        {
        mat2x2  = matgen2[ , c(i,j) ]

        k       = PAIRINDEX_plus( i, j, n )

        b   = qpoint - centerrot[k,1:2]

        alpha   = as.double( solve( mat2x2, b ) )

        absalpha    = abs(alpha)

        #cat( "===============", "  iter ", iter, "  =================\n" )
        #cat( "i =", i, "  j =", j,  "   alpha =", alpha, "   center=",  centerrot[k,1:2],   "  d2 =", sum(b^2), "  z_delta =", centerrot[k,3]-baserot[3], '\n' )
        #cat( "    det =", det2x2(mat2x2), '\n' )

        if( all( absalpha <= 1/2 ) )
            {
            #   qpoint is inside the pgram !   We got it !
            success = TRUE
            break
            }

        #   move to the next pgram2
        if( which.max(absalpha) == 1 )
            {
            #   change i
            if( 0 < alpha[1] )
                {
                i   = idx_dec[i]                    # decrement i
                if( i == j )    j   = idx_dec[j]    # decrement j too !
                }
            else
                {
                i   = idx_inc[i]                    # increment i
                if( i == j )    j   = idx_inc[j]    # increment j too !
                }
            }
        else
            {
            #   change j
            if( 0 < alpha[2] )
                {
                j   = idx_inc[j]                    # increment j
                if( j == i )    i   = idx_inc[i]    # increment i too !
                }
            else
                {
                j   = idx_dec[j]                    # decrement j
                if( j == i )    i   = idx_dec[i]    # decrement i too !
                }
            }
        }

    if( ! success )
        {
        log_level( ERROR, "Reached maximum iterations: %d.", maxiters )
        return(NULL)
        }


    out = list()
    out$idxpair     = c(i,j)
    out$alpha       = alpha + 1/2   # from centered to uncentered
    out$iters       = iter

    #print( out )

    #cat( "findpgram2D(). timesearch =", gettime() - time_start, '\n' )

    return( out )
    }



if( FALSE )
{

#   given a tiling by pgrams of a part of the plane, and a query point
#   find the pgram that contains the query point
#
#   centerrot   N*(N-1) x 3 matrix with pgram centers,
#               only the first 2 columns are used, the 3rd z-coordinate might be used in future testing
#   baserot     3-vector with query point, 3rd z-coord not used
#               the goal is to find the pgram that contains this point
#   matgen2     2 x N matrix of generators
#   idxpair0    integer pair to start the search

findpgram2D_v2   <- function( centerrot, baserot, matgen2, idxpair0, tol=5.e-9 )
    {
    time_start = gettime()

    n   = ncol(matgen2)

    if( nrow(centerrot) != n*(n-1) )
        {
        log_level( ERROR, "nrow(centerrot) = %g != %g = n*(n-1).  n=%g", nrow(centerrot), n*(n-1), n )
        return(NULL)
        }

    qpoint  = baserot[1:2]

    #   at each point in the iteration, the state is determined by 2 things:
    #       the current pgram2, given by i and j
    #       ppoint, which is a point inside the current pgram, usually a boundary point but initially the center
    #           when a boundary point, it is shared with the previous pgram

    #   use variable i and j to abbreviate 2 ints in idxpair0
    i   = idxpair0[1]
    j   = idxpair0[2]

    k       = PAIRINDEX_plus( i, j, n )
    ppoint  = centerrot[k,1:2]

    maxiters    = as.integer( max( 0.5*n, 50 ) )
    success     = FALSE

    #   make increment and decrement vectors
    idx_inc = c(2:n,1)
    idx_dec = c(n,1:(n-1))

    for( iter in 1:maxiters )
        {
        mat2x2  = matgen2[ , c(i,j) ]
        mat2x2_inv  = solve( mat2x2 )

        k   = PAIRINDEX_plus( i, j, n )

        # test whether qpoint is inside pgram (i,j)
        vec = qpoint - centerrot[k,1:2]

        alpha   = as.double( mat2x2_inv %*% vec )       # solve( mat2x2, vec )

        cat( "===============", "  iter ", iter, "  =================\n" )
        dist    = sqrt( sum((qpoint - ppoint)^2) )
        cat( "i =", i, "  j =", j,  "   alpha =", alpha, "   center=",  centerrot[k,1:2],   "  dist =", dist, "  unitize(qpoint-ppoint) =", (qpoint-ppoint)/dist, '\n' )

        absalpha    = abs(alpha)

        if( all( absalpha <= 1/2 ) )
            {
            #   qpoint is inside the pgram !   We got it !
            success = TRUE
            break
            }

        #   move to the next pgram2

        #   transform point in pgram and direction (qpoint - ppoint) to the unit square
        b   = mat2x2_inv %*% (ppoint - centerrot[k,1:2])    # b is in square [-1/2,1/2]^2
        v   = mat2x2_inv %*% (qpoint - ppoint)              # v points into the interior of the square

        #   snap to +-1/2 if within tol
        snap    = abs(abs(b) - 1/2) <= tol
        b[snap] = round( b[snap] + 1/2 ) - 1/2

        cat( "b before =", b, "  v =", v,  "   ppoint =", ppoint, '\n' )


        b   = squareadvance( b, v ) # new b is on boundary of square

        ppoint  = as.double( mat2x2 %*% b )  +  centerrot[k,1:2]

        cat( "b after =", b, "   ppoint =", ppoint, '\n' )

        #return(NULL)

        if( abs(b[1]) == 0.5 )
            {
            #   left or right edge, change i
            if( b[1] == 0.5 )
                {
                i   = idx_dec[i]                    # decrement i
                if( i == j )    j   = idx_dec[j]    # decrement j too !
                }
            else
                {
                i   = idx_inc[i]                    # increment i
                if( i == j )    j   = idx_inc[j]    # increment j too !
                }
            }
        else
            {
            #   bottom or top edge, change j
            if( b[2] == 0.5 )
                {
                j   = idx_inc[j]                    # increment j
                if( j == i )    i   = idx_inc[i]    # increment i too !
                }
            else
                {
                j   = idx_dec[j]                    # decrement j
                if( j == i )    i   = idx_dec[i]    # decrement i too !
                }
            }
        }

    if( ! success )
        {
        log_level( ERROR, "Reached maximum iterations: %d.", maxiters )
        return(NULL)
        }


    out = list()
    out$idxpair     = c(i,j)
    out$alpha       = alpha + 1/2   # from centered to uncentered
    out$iters       = iter

    timesearch = gettime() - time_start

    cat( "findpgram2D_v2(). timesearch =", timesearch, '\n' )

    return( out )
    }
}



#   in this version, use brute-force search to find the containing pgram
#
#   centerrot       N*(N-1) x 3 matrix with pgram centers,
#                   the 3rd z-coordinate is used to skip pgrams that are too far "below" the basepoint
#   baserot         3-vector with query point, 3rd z-coord is used too
#                   the goal is to find the pgram that contains this point
#   idxpair_plus    N*(N-1) x 2 integer matrix of 1-based generator indexes, including antipodal data
#   genrot          3 x N matrix of rotated generators

findpgram2D_v3  <- function( centerrot, baserot, idxpair_plus, genrot )
    {
    #print( centerrot )
    #cat( "baserot =", baserot, '\n' )

    res = .Call( C_findpgram2D, centerrot, baserot, idxpair_plus, genrot )

    if( res[[1]] < 0 )  return(NULL)

    k       = res[[1]] + 1L     # 0-based to 1-based
    alpha   = res[[2]]

    out = list()
    out$idx         = k
    out$idxpair     = idxpair_plus[k, ]
    out$alpha       = alpha + 1/2           # from centered to uncentered

    #   out$iters       = k

    return( out )
    }





#   base            base of ray, uncentered
#   dir             direction of ray
#   center          center of of symmetry
#   matgen          3 x M matrix of generators
#   pgramdf         data frame with idxpair, center, and other variables.   not extended with antipodal data

#   returns a list with these items:
#       idxpair     indexes of the generators of the pgram, or both NA if the ray intersects a pole
#       alpha       2 coords in [0,1] for point inside the pgram
#       iters       # of iterations, or 0 when the ray intersect a pgram that intersects a pole
#       tmax        parameter of ray when it intersects the surface, or NA in case of failure
#       point       intersection of ray with the surface, or NA in case of failure

findpgram3D   <- function( base, dir, center, matgen, pgramdf )
    {
    #   extend with antipodal data
    pgramdf = pgramdf_plus( pgramdf )

    base_centered = base - center

    success     = FALSE

    for( k in 1:nrow(pgramdf) )
        {
        idxpair     = pgramdf$idxpair[k, ]

        mat = cbind( matgen[ , idxpair ], -dir )

        y   = solve( mat, base_centered - pgramdf$center[k, ] ) #;        print(y)

        alpha   = y[1:2]
        tau     = y[3]

        if( all( abs(alpha) <= 0.5 )  &&  0 < tau )
            {
            success = TRUE
            break
            }
        }

    if( ! success )     return(NULL)

    out = list()
    out$idx     = k
    out$idxpair = idxpair
    out$alpha   = alpha + 1/2       # centered to uncentered

    #   print( out )

    return( out )
    }











#   this is for use with pgramdf_plus()
#   this is only valid if 1 <= i , j <= n.  But if i==j it returns NA_integer_
PAIRINDEX_plus  <- function( i, j, n )
    {
    if( i < j )
        out = (i-1L)*n - ((i)*(i+1L)) %/% 2L + j
    else if( j < i )
        out = (j-1L)*n - ((j)*(j+1L)) %/% 2L + i  + (n*(n-1L)) %/% 2L
    else
        out = NA_integer_

    return( out )
    }



if( FALSE )
{
#   b   point in the square [-1/2,1/2]^2
#       usually it is on the boundary, or the center c(0,0) - not verified
#   v   non-zero vector pointing into the interior of the square  -  not verified
#   tol tolerance for snapping output to the boundary
#
#   returns point where the ray     b + t*v     intersects the boundary
#   returns NULL in case of problem

squareadvance   <- function( b, v, tol=5.e-8 )
    {
    #   check b
    #absb    = abs(b)
    #ok  =  1 <= sum( absb==1/2 )  &&  all( absb <= 1/2 )
    #if( ! ok )  return(NULL)

    tvec = c( (-0.5 - b)/v, (0.5 - b)/v )

    #   change any non-positive or NaN entries to Inf
    tvec[ tvec <= 0 | is.nan(tvec) ] = Inf

    #   find the minimum entry
    tmin    = min( tvec )
    if( ! is.finite(tmin) ) return(NULL)

    out = b + tmin*v

    #   snap to +-1/2 if within tol
    snap        = abs(abs(out) - 1/2) <= tol
    out[snap]   = round( out[snap] + 1/2 ) - 1/2

    #   check output
    #ok  = all( abs(out) <= 1/2 )
    #if( ! ok )  return(NULL)

    cat( "squareadvance(). tmin =", tmin, "   out =", out, '\n' )

    return( out )
    }
}




#   theta       query angle
#   thetaend    2 angles forming the endpoints of an arc.  The shorter arc is the one intended.
#
#   returns TRUE if theta is on the arc

on_arc  <- function( theta, thetaend )
    {
    #   rotate both endpoints to put theta at 0,
    #   and both endpoints in [-pi,pi)
    thetaend    = ( (thetaend - theta + pi) %% (2*pi) ) - pi

    #   if( any(thetaend == 0) )    return(TRUE)    # theta is an endpoint

    if( pi <= abs(thetaend[2] - thetaend[1]) )
        #   0 is in the wrong arc
        return(FALSE)

    #   check whether the endpoints are on opposite sides of 0
    return( prod(thetaend) <= 0 )
    }




#   mask    a logical mask
#
#   returns TRUE  iff  the set of TRUE values in mask[] is contiguous in a cyclic sense

is_contiguousmask    <- function( mask )
    {
    subs    = which( mask )

    if( all( diff(subs) == 1L ) )
        return(TRUE)

    #   try the complements
    subs    = which( ! mask )

    return( all( diff(subs) == 1L ) )
    }


det2x2  <- function( mat2x2 )
    {
    return( mat2x2[1,1]*mat2x2[2,2] - mat2x2[1,2]*mat2x2[2,1] )
    }



if( FALSE )
{
#   idxpair     count x 2 integer matrix of pgrams to plot

plotpgrams2D <- function( x, base, direction, idxpair, winrad=c(0.001,0.001) )
    {
    center  = getcenter(x)

    base_centered = base - center

    #   from current direction, compute a 2x3 projection matrix
    frame3x3    = goodframe3x3( direction )

    ok  = is.matrix(idxpair)  &&  is.integer(idxpair)  &&  ncol(idxpair)==2

    pgramdf = allpgramcenters2trans(x)
    if( is.null(pgramdf) )  return(NULL)

    matsimple   = getsimplified( x$matroid )

    matgen  = getmatrix( matsimple )

    m       = ncol(matgen)

    matproj     = t( frame3x3[ , 1:2 ] )

    matgen2 = matproj %*% matgen

    temp    = pgramdf$center %*% frame3x3

    centerrot   = .Call( C_extend_antipodal, temp )

    baserot     = as.double( base_centered %*% frame3x3 )
    cat( "baserot = ", baserot, '\n' )

    pgrams  = nrow(idxpair)

    #   allocate 3D array for all the vertices
    vertex  = array( 0, c(pgrams,4,2) )

    for( i in 1:pgrams )
        {
        idx     = idxpair[i, ]

        k   = PAIRINDEX_plus( idx[1], idx[2], m )

        center  = centerrot[k,1:2]

        #   get the 2 generators of the pgram
        gen    = matgen2[  , idx ]

        vertex[i,1, ]   = center - 0.5*gen[ ,1] - 0.5*gen[ ,2]
        vertex[i,2, ]   = center - 0.5*gen[ ,1] + 0.5*gen[ ,2]
        vertex[i,3, ]   = center + 0.5*gen[ ,1] + 0.5*gen[ ,2]
        vertex[i,4, ]   = center + 0.5*gen[ ,1] - 0.5*gen[ ,2]
        }


    if( FALSE )
        {
        #   pgrams are typically very thin, so apply whitening
        vertex_xy   = cbind( as.double(vertex[ , ,1]), as.double(vertex[ , ,2]) ) ;
        print( vertex_xy )

        #   center vertex_xy
        vertex_xy   = vertex_xy - matrix( colMeans(vertex_xy), nrow(vertex_xy), ncol(vertex_xy), byrow=TRUE )
        print( vertex_xy )


        res = base::svd( vertex_xy, nu=0, nv=2 )

        if( ! all( 0 < res$d ) )
            {
            log_level( ERROR, "Internal error. vertex matrix is invalid." )
            return(NULL)
            }

        #   method == "ZCA" )
        #   calculate the "whitening", or "sphering" matrix, which is MxM
        W   = res$v %*% diag( 1/res$d ) %*% t(res$v)        ; print( W )

        vertex  = vertex_xy %*% t(W)

        dim(vertex)     = c(pgrams,4,2)
        }


    xlim    = range( vertex[ , ,1] )
    ylim    = range( vertex[ , ,2] )

    xlim = baserot[1] + c(-1,1) * winrad[1]
    ylim = baserot[2] + c(-1,1) * winrad[2]


    plot( xlim, ylim, type='n', las=1, asp=1, lab=c(10,10,7), xlab='x', ylab='y' )
    grid( lty=1 )
    abline( h=0, v=0 )

    for( i in 1:pgrams )
        {
        col = ifelse( i==pgrams, 'lightyellow', NA )

        polygon( vertex[i, ,1], vertex[i, ,2], col=col )

        idx     = idxpair[i, ]
        k   = PAIRINDEX_plus( idx[1], idx[2], m )
        points( centerrot[k,1], centerrot[k,2], pch=20 )
        }

    points( baserot[1], baserot[2] )


    return( invisible(TRUE) )
    }
}




#   section() compute intersection of 2-transition sphere and plane(s)
#
#   x           a zonohedron object
#   normal      a non-zero numeric vector of length 3, the normal of all the planes
#   beta        a vector of plane-constants.  The equation of plane k is: <x,normal> = beta[k]
#
#   value   a list of data.frames with length = length(beta).

section2trans <- function( x, normal, beta, invert=FALSE, plot=FALSE, tol=1.e-12, ... )
    {
    if( ! inherits(x,"zonohedron") )
        {
        log_level( ERROR, "Argument x is invalid.  It's not a zonohedron." )
        return(NULL)
        }


    ok  = is.numeric(normal)  &&  length(normal)==3  &&  all( is.finite(normal) )  &&  any( normal!=0 )
    if( ! ok )
        {
        log_level( ERROR, "normal is invalid. It must be a non-zero numeric vector of length 3, and all entries finite." )
        return(NULL)
        }

    ok  = is.numeric(beta)  &&  0<length(beta)   &&  all( is.finite(beta) )
    if( ! ok )
        {
        log_level( ERROR, "beta is invalid. It must be a numeric vector of positive length, and all entries finite." )
        return(NULL)
        }

    cnames  = names(normal) # save these
    dim(normal) = NULL

    matsimple   = getsimplified(x$matroid)
    matgen      = getmatrix( matsimple )
    gndgen      = getground( matsimple )

    numgen  = ncol( matgen )    # so matgen is 3 x numgen


    #   make increment and decrement vectors
    ij_inc  = c(2:numgen,1L)
    ij_dec  = c(numgen,1:(numgen-1))



    #   for each generator, compute radius of the projection of the generator onto the line generated by normal
    normalgen   = as.numeric( normal %*% matgen )
    radiusgen   = 0.5 * abs(normalgen)              # so length(radiusgen) = numgen

    pgramdf = allpgramcenters2trans(x) #;    print( str(pgramdf) )
    if( is.null(pgramdf) )  return(NULL)

    #   the radius of a pgram is simply the sum of the radii of the generators
    myfun   <- function( pair )    { sum( radiusgen[ pair ] ) }

    radiuspgram = apply( pgramdf$idxpair, 1L, myfun ) #;    print( str(radiuspgram) )


    #   dot products of the pgram centers and the normal vector
    cn      = as.numeric( pgramdf$center %*% normal )

    #   append both center and cn with the antipodal data
    center  = rbind( pgramdf$center, -pgramdf$center )
    cn      = c( cn, -cn )
    idxpair = rbind( pgramdf$idxpair,  pgramdf$idxpair[ ,2:1] )

    numpgrams   = length( cn )  # includes antipodal data

    #   compute range of normal over each pgram
    cnneg   = cn - radiuspgram  # the recycling rule is used here
    cnpos   = cn + radiuspgram  # the recycling rule is used here

    #   translate beta to the centered zonogon
    cent.norm = sum( x$center * normal )
    betacentered   = as.numeric(beta) - cent.norm  #; print(betacentered)


    #   load the coefficients that traverse a parallelogram first by generator 1 and then generator 2.
    #   we call this the "standard" order
    vertexcoeff = matrix( c( -0.5,-0.5, 0.5,-0.5, 0.5,0.5, -0.5,0.5), 4, 2, byrow=TRUE )

    #   scratch vector that holds dot product of the parallelogram vertices with the normal
    value   = numeric( 4 )

    next4   = c( 2:4, 1L )

    #   the intersection of the plane and the polyhedral surface is a polygon
    #   make matrix to hold all vertices of the polygons
    vertex      = matrix( NA_real_, nrow=nrow(center), ncol=3 )
    adjacent    = integer( numpgrams )


    out         = vector( length(beta), mode='list' )
    names(out)  = sprintf( "normal=%g,%g,%g. beta=%g", normal[1], normal[2], normal[3], beta )

    for( k in 1:length(beta) )
        {
        beta_k  = betacentered[k]

        #   find indexes of all pgrams whose *interiors* intersect this plane
        #   NOTE: If a pgram intersects the plane only in a vertex or edge, it will not be found !
        #   perhaps fix this later ?
        maskinter   = cnneg < beta_k  &  beta_k < cnpos

        indexvec    = which( maskinter )  #; print( indexvec )

        if( length(indexvec) < 3 )
            {
            #   plane does not intersect the zonohedron, there is no section
            #   or the section is a degenerate single point or an edge
            # out[[k]]    = list( beta=beta[k],  section=matrix( 0, 0, 3 ) )
            out[[k]]        = data.frame( row.names=integer(0) )
            out[[k]]$point  = matrix(0,0,3)
            next
            }

        vertex[ , ] = NA_real_      # clear vertex
        adjacent[ ] = NA_integer_   # clear adjacent

        if( invert )    pcube = matrix( NA_real_, length(indexvec), numgen )


        #   length(indexvec) is the # of pgrams that the plane intersects
        for( ii in 1:length(indexvec) )
            {
            #   idx is a row index into center[] and idxpair and vertex[]
            idx     = indexvec[ii]

            #   find point where the plane intersects an edge of the pgram
            pair    = idxpair[ idx, ]   #sort( idxpair[ idx, ] )
            i       = pair[1]
            j       = pair[2]

            #   compute the values at the vertices in the 'standard' order
            #for( ii in 1:4 )
            #    value[ii] = vertexcoeff[ii,1] * normalgen[i]  +  vertexcoeff[ii,2] * normalgen[j]  +  cn[idx] - beta_k

            value   = as.double( vertexcoeff %*% normalgen[pair] )  +  cn[idx] - beta_k

            #   traverse value[] and check that there is exactly 1 transition from neg to non-neg
            #itrans  = integer(0)
            #for( ii in 1:4 )
            #    {
            #    if( value[ii]<0  &&  0<=value[ next4[ii] ] )
            #        itrans = c(itrans,ii)   #;  count = count+1L
            #    }

            itrans  = which( value<0  &  0<=value[next4] )

            if( length(itrans) != 1 )
                {
                log_level( ERROR, "Internal Error. pgram idx=%d.  value[] has %d transitions, but expected 1.", idx, length(itrans) )
                next
                }

            #   we have found the right edge of the pgram
            #   find the intersection of the edge and the plane
            i1 = itrans
            i2 = next4[i1]

            v1  = as.numeric( matgen[ ,pair] %*% vertexcoeff[i1, ] )
            v2  = as.numeric( matgen[ ,pair] %*% vertexcoeff[i2, ] )

            lambda = value[i2] / (value[i2] - value[i1])

            vertex[idx, ]   = center[idx, ] + lambda*v1  +  (1-lambda)*v2


            #   find the pgram on the other side of this edge
            #   this adjacent pgram must also be in indexvec
            if( itrans == 1 )
                {
                j   = ij_dec[j]                     # decrement j
                if( j == i )    i   = ij_dec[i]     # decrement i too !
                alpha   = c(1-lambda,0)
                }
            else if( itrans == 2 )
                {
                i   = ij_dec[i]                     # decrement i
                if( i == j )    j   = ij_dec[j]     # decrement j too !
                alpha   = c(1,1-lambda)
                }
            else if( itrans == 3 )
                {
                j   = ij_inc[j]                     # increment j
                if( j == i )    i   = ij_inc[i]     # increment i too !
                alpha   = c(lambda,1)
                }
            else if( itrans == 4 )
                {
                i   = ij_inc[i]                     # increment i
                if( i == j )    j   = ij_inc[j]     # increment j too !
                alpha   = c(0,lambda)
                }

            adjacent[idx]   = PAIRINDEX_plus( i, j, numgen )

            if( invert )    pcube[ii, ]     = pcubefromdata( pair, alpha, numgen )

            #cat( "maskinter[ adjacent[idx] ] =", maskinter[ adjacent[idx] ], "    value[i2] =", value[i2], '\n' )

            if( abs(value[i2])<=tol  &&  ! maskinter[ adjacent[idx] ] )
                {
                #   degenerate case, try to fix it
                #cat( "fixing itrans =", itrans, '\n' )
                if( itrans == 1 )
                    {
                    if( i == pair[1] )    i   = ij_dec[i]     # decrement i too !
                    }
                else if( itrans == 3 )
                    {
                    if( i == pair[1] )    i   = ij_inc[i]     # increment i too !
                    }

                # this still might be an invalid pgram; if so it will be caught later
                adjacent[idx]   = PAIRINDEX_plus( i, j, numgen )
                }
            }

        #print( indexvec )
        #print( adjacent )


        #   put the non-trivial rows of vertex[,] in proper order,
        #   using the adjacent[] index vector

        success = TRUE

        indexordered    = rep( NA_integer_, length(indexvec) )

        #   since the *interiors* of all these pgrams intersect the plane,
        #   we can start this iteration at any one of them
        indexordered[1] = indexvec[1]
        for( i in 2:length(indexordered) )
            {
            indexordered[i] = adjacent[ indexordered[i-1] ]

            if( is.na(indexordered[i]) )
                {
                log_level( WARN, "Internal Error.  bad adjacent logic.  i=%d.", i )
                success = FALSE
                break
                }

            if( indexordered[i] == indexordered[1] )
                {
                #   back to the starting pgram, premature stop, back off !
                #   indexordered[i] = NA_integer_
                i = i - 1L
                break
                }
            }

        if( FALSE )
            {
            print( indexvec )
            print( vertex )
            print( adjacent )
            print( indexordered )
            }

        if( i < length(indexordered) )
            {
            log_level( WARN, "The plane for beta=%g intersects the surface in more than 1 polygon.", beta[k] )    #; but only 1 polygon is being returned
            #log.string( ERROR, "The plane for beta=%g intersects the surface in more than 1 polygon.", beta[k] )
            success = FALSE
            #indexordered    = indexordered[1:i]     #   trim away the excess
            }


        if( ! success )
            {
            #   assign the empty section
            #   out[[k]]    = list( beta=beta[k],  section=matrix( 0, 0, 3 ) )
            out[[k]]        = data.frame( row.names=integer(0) )
            out[[k]]$point  = matrix(0,0,3)
            next
            }



        #   extract the polygon
        poly    = vertex[indexordered, ]

        #   the poly coordinates are centered, so add back the zono center
        res = .Call( C_plusEqual, poly, x$center, 1L )   # changes poly in place
        if( is.null(res) )  return(NULL)

        df  = data.frame( row.names=1:nrow(poly) )

        df$point    = poly

        gndpair     = gndgen[ idxpair[ indexordered, ] ]
        dim(gndpair)    = c( length(indexordered), 2 )
        df$gndpair  = gndpair

        #gndpairadj      = gndgen[ idxpair[ adjacent[indexordered], ] ]
        #dim(gndpairadj) = c( length(indexordered), 2 )
        #df$gndpairadj   = gndpairadj

        if( invert )
            {
            perm        = order( indexordered )
            perm[perm]  = 1:length(perm)  # invert perm
            df$pcube    = invertcubepoints( x, pcube[ perm, ] )

            #print( indexordered )
            #print( perm )
            }

        out[[k]]    = df

        if( FALSE  &&  invert )
            {
            #   test the inversion
            #print( df$pcube )

            matorig = getmatrix( x$matroid )

            delta   = df$pcube %*% t(matorig)  -  df$point

            # cat( "range(delta)=", range(delta), '\n' )

            delta   = rowSums( abs(delta) )

            #print( t(matorig) )

            if( any( tol < delta, na.rm=TRUE ) )
                {
                #print( delta )
                #log.string( WARN, "Inversion test failed.  max(delta)=%g > %g=tol",
                #                    max(delta,na.rm=TRUE), tol )
                }

            out[[k]]$delta   = delta
            }
        }




    if( plot )
        {
        if( ! requireNamespace( 'rgl', quietly=TRUE ) )
            {
            log_level( WARN, "Plotting cannot be done, because package 'rgl' is not installed. " )
            }
        else if( rgl::cur3d() == 0 )
            {
            log_level( WARN, "Cannot add section to plot, because there is no rgl window open." )
            }
        else
            {
            for( k in 1:length(beta) )
                {
                xyz = out[[k]]$point
                rgl::polygon3d( xyz[ ,1], xyz[ ,2], xyz[ ,3], fill=FALSE, col='red' )
                #   rgl::points3d( xyz[ ,1], xyz[ ,2], xyz[ ,3], col='red', size=13, point_antialias=TRUE )
                }
            }
        }

    return( invisible(out) )
    }





transitionsdf <- function( x, trans2=TRUE )
    {
    if( ! inherits(x,"zonohedron") )
        {
        log_level( ERROR, "Argument x is invalid.  It is not a zonohedron." )
        return(NULL)
        }

    metrics   = getmetrics2trans( x, angles=FALSE )
    if( is.null(metrics) )   return(NULL)

    ok  = is.finite(metrics$starshaped)  &&  metrics$starshaped
    if( ! ok )
        {
        log_level( ERROR, "The zonohedron x is invalid.  The 2-transition surface is NOT strictly starshaped." )
        return(NULL)
        }

    #print( str( metrics$pgramdf ) )

    deficient   = metrics$pgramdf$deficient

    #deficientcount  = sum(deficient)

    #dfacets = 2 * deficientcount
    #cat( "deficient parallelograms:         ", dfacets, "   [fraction=", dfacets/facets, ']\n' )
    #cat( "deficient area:                   ", metrics$areadeficient, "   [fraction=", metrics$areadeficient/res$area, ']\n' )
    #cat( "deficient volume:                 ", metrics$volumedeficient, "   [fraction=", metrics$volumedeficient/res$volume, ']\n' )

    #thickness = metrics$volumedeficient / metrics$areadeficient
    #cat( "deficient thickness (mean):       ", thickness, '\n' )

    if( any(deficient) )
        {
        #   compute transitions, but only for the deficient parallelograms
        pgramdef    = metrics$pgramdf[ deficient, ]

        dat         = boundarypgramdata( x, pgramdef$gndpair  )
        
        #   dat now has these columns: gndpair, hyperplaneidx, center, transitions
        #   dat$gndpair = pgramdef$gndpair
        
        #   now add more columns
        dat$area    = pgramdef$area
        dat$deficit = pgramdef$deficit

        #   remove unneeded columns
        dat$hyperplaneidx   = NULL
        dat$center          = NULL
        }
    else
        dat = NULL

    #print( str(dat) )

    notdef  = ! deficient

    if( trans2  &&  any(notdef) )
        {
        #   append more rows to dat
        #   append the rows to dat that are *not* deficient and therefore have 2 transitions
        dat2    = data.frame( row.names=1:sum(notdef) )
        dat2$gndpair        = metrics$pgramdf$gndpair[ notdef, ]
        dat2$transitions    = 2L
        dat2$area           = metrics$pgramdf$area[ notdef ]
        dat2$deficit        = metrics$pgramdf$deficit[ notdef ]

        #print( str(dat2) )

        dat = rbind( dat, dat2 )
        }

    if( is.null(dat) )
        {
        log_level( ERROR, "No rows to return." )
        return(NULL)
        }


    #   breakdown the deficient parallelograms by transitions
    tunique = sort( unique( dat$transitions )  )    #; print( tunique )
    m   = length(tunique)

    pgcount = integer(m)

    area.r              = matrix( 0, nrow=m, ncol=2 )
    colnames(area.r)    = c( 'min', 'max' )

    deficit.r           = matrix( 0, nrow=m, ncol=2 )
    colnames(deficit.r) = c( 'min', 'max' )

    area.s              = numeric( m )

    example = character(m)

    for( k in 1:m )
        {
        datsub      = dat[ dat$transitions == tunique[k], ]
        pgcount[k]  = 2L*nrow( datsub )
        area.r[k, ] = range( datsub$area )
        area.s[k]   = 2*sum( datsub$area )

        deficit.r[k, ]  = range( datsub$deficit )

        i  = which.max( datsub$area )
        gndpair     = datsub$gndpair[i, ]
        example[k]  = sprintf( "{%d,%d}", gndpair[1], gndpair[2] )
        }

    out = data.frame( row.names = c(1:m,"Totals") )

    out$transitions     = c(tunique,NA)
    out$parallelograms  = c( pgcount, sum(pgcount) )
    out$area            = rbind( area.r, c(NA,NA) )
    out$area.sum        = c(area.s,sum(area.s))
    out$deficit         = rbind( deficit.r, c(NA,NA) )
    out$example         = c( example, '' )

    return( out )
    }




plothighertrans <- function( x, abalpha=1, defcol='green', defalpha=0, ecol=NA,
                                connections=FALSE, bgcol="gray40", both=TRUE, ... )
    {
    if( ! inherits(x,"zonohedron") )
        {
        log_level( ERROR, "Argument x is invalid.  It's not a zonohedron." )
        return(NULL)
        }


    if( ! requireNamespace( 'rgl', quietly=TRUE ) )
        {
        log_level( ERROR, "Package 'rgl' cannot be loaded. It is required for plotting the zonohedron." )
        return(FALSE)
        }

    if( abalpha<=0  &&  defalpha<=0  )
        {
        log_level( WARN, "abalpha=%g and defalpha=%g is invalid.", abalpha, defalpha )
        return( FALSE )
        }

    metrics = getmetrics2trans( x, angles=FALSE, tol=1.e-12 )

    if( is.null(metrics) )
        return(FALSE)

    ok  = is.finite(metrics$starshaped)  &&  metrics$starshaped
    if( ! ok )
        {
        log_level( ERROR, "The zonohedron x is invalid.  The 2-transition surface is NOT strictly starshaped." )
        return(NULL)
        }

    facesok  = ! is.null(metrics$pgramdf$deficient)     #; cat( 'facesok ', facesok, '\n' )

    if( ! facesok )
        {
        log_level( WARN, "Cannot draw faces because the parallelogram data are not available." )
        return( FALSE )
        }


    center  = x$center
    white   = 2 * center


    #   start 3D drawing
    rgl::bg3d( color=bgcol )

    #cube    = rgl::scale3d( rgl::cube3d(col="white"), center[1], center[2], center[3] )
    #cube    = rgl::translate3d( cube, center[1], center[2], center[3]  )

    rgl::points3d( 0, 0, 0, col='black', size=10, point_antialias=TRUE )
    rgl::points3d( white[1], white[2], white[3], col='white', size=10, point_antialias=TRUE )
    rgl::points3d( center[1], center[2], center[3], col='gray50', size=10, point_antialias=TRUE )

    #   exact diagonal
    rgl::lines3d( c(0,white[1]), c(0,white[2]), c(0,white[3]), col=c('black','white'), lwd=3, lit=FALSE )


    matsimp = getsimplified( x$matroid )

    matgen  = getmatrix(matsimp)
    numgen  = ncol(matgen)

    #gndgen  = getground(matsimp)



    pgramdf = metrics$pgramdf

    deficient   =  pgramdf$deficient

    deficientcount  = sum( deficient )

    if( deficientcount == 0 )
        {
        log_level( WARN, "Cannot draw because none of the 2-transition facets are deficient." )
        return( FALSE )
        }

    drawedges   = ! is.na(ecol)

    #   extract only the subset of deficient rows
    pgramdf = pgramdf[ deficient, ]  #; print( pgramdf )

    # pgramdf$alphavec    = alphavec

    step    =  4
    quadmat = matrix( 0, nrow=step*deficientcount, ncol=3 )


    if( 0 < abalpha )
        {
        #   draw filled abundant pgrams

        #   use the same colorkey as Scott Burns
        colorkey    = c( "black", "white", "black", "#8c0000", "black", "#ffff19", "black", "#0082c8", "black", "#dcbeff", "green" )

        #   get bounddf in order to get the # of transitions
        bounddf = boundarypgramdata( x, pgramdf$gndpair )

        colvec      = character( step*deficientcount )
        #alphavec    = numeric( step*deficientcount )

        for( i in 1:deficientcount )
            {
            center  = pgramdf$centermax[i, ]

            edge    = matgen[  , pgramdf$idxpair[i, ] ]  # 3x2 matrix

            k       = step*(i-1)

            quadmat[k+1, ] = center - 0.5 * edge[ , 1] - 0.5*edge[ , 2]
            quadmat[k+2, ] = center - 0.5 * edge[ , 1] + 0.5*edge[ , 2]
            quadmat[k+3, ] = center + 0.5 * edge[ , 1] + 0.5*edge[ , 2]
            quadmat[k+4, ] = center + 0.5 * edge[ , 1] - 0.5*edge[ , 2]

            tcount  = min( bounddf$transitions[i], length(colorkey) )

            colvec[ (k+1):(k+4) ]   = colorkey[ tcount ]       # pgramdf$colvec[i]

            #alphavec[ (k+1):(k+4) ] = pgramdf$alphavec[i]
            }

        xyz = .Call( C_sumMatVec, quadmat, x$center, 1L )
        rgl::quads3d( xyz, col=colvec, alpha=abalpha, lit=FALSE  )              # quad filled

        if( drawedges )
            rgl::quads3d( xyz, col=ecol, lwd=2, front='lines', back='lines', lit=FALSE  )   # quad edges

        if( both )
            {
            xyz = .Call( C_sumMatVec, -quadmat, x$center, 1L )
            rgl::quads3d( xyz, col=colvec, alpha=abalpha, lit=FALSE )

            if( drawedges )
                rgl::quads3d( xyz, col=ecol, lwd=2, front='lines', back='lines', lit=FALSE  )   # quad edges
            }

        #   rgl::legend3d( "top", c("2D Points", "3D Points"), cex=0.75,  pch = c(1, 16)  )

        #   rgl::title3d('main', 'sub', 'xlab', 'ylab', 'zlab')
        }


    if( 0 < defalpha )
        {
        #   draw filled deficient pgrams
        for( i in 1:deficientcount )
            {
            center  = pgramdf$center[i, ]

            edge    = matgen[  , pgramdf$idxpair[i, ] ]  # 3x2 matrix

            k       = step*(i-1)

            quadmat[k+1, ] = center - 0.5 * edge[ , 1] - 0.5*edge[ , 2]
            quadmat[k+2, ] = center - 0.5 * edge[ , 1] + 0.5*edge[ , 2]
            quadmat[k+3, ] = center + 0.5 * edge[ , 1] + 0.5*edge[ , 2]
            quadmat[k+4, ] = center + 0.5 * edge[ , 1] - 0.5*edge[ , 2]
            }

        xyz = .Call( C_sumMatVec, quadmat, x$center, 1L )
        rgl::quads3d( xyz, col=defcol, alpha=defalpha, lit=FALSE  )              # quad filled

        if( drawedges )
            rgl::quads3d( xyz, col=ecol, lwd=2, front='lines', back='lines', lit=FALSE  )   # quad edges

        if( both )
            {
            xyz = .Call( C_sumMatVec, -quadmat, x$center, 1L )
            rgl::quads3d( xyz, col=defcol, alpha=defalpha, lit=FALSE )

            if( drawedges )
                rgl::quads3d( xyz, col=ecol, lwd=2, front='lines', back='lines', lit=FALSE  )   # quad edges
            }
        }




    #   for each pair of pgrams in the zonohedron and the 2-transition polyhedron,
    #   draw a segment connecting their 2 centers
    #   these are matching deficient and abundant parallelograms

    if( connections )
        {
        #   draw segments between the deficient pgrams in the 2-transition surface,
        #   and the corresponding pgrams in the zonohedron boundary
        #pgramdf = metrics$pgramdf
        #centerdef   = pgramdf$center[deficient, ]   #metrics$signsurf *     metrics$
        #centermax   = pgramdf$centermax[ deficient, ]

        centermat   = rbind( pgramdf$center, pgramdf$centermax )

        mat = matrix( 1:nrow(centermat), nrow=deficientcount, ncol=2 )
        idx = as.integer( t(mat) )
        centermat   = centermat[ idx, ]     #; print( centermat )

        xyz = .Call( C_sumMatVec, centermat, x$center, 1L )
        rgl::segments3d( xyz, col='black'  )
        rgl::points3d( xyz, col=c('yellow','black'), size=5, point_antialias=TRUE  )

        if( both )     # both
            {
            xyz = .Call( C_sumMatVec, -centermat, x$center, 1L )
            rgl::segments3d( xyz, col='black' )
            rgl::points3d( xyz, col=c('yellow','black'), size=5, point_antialias=TRUE  )
            }
        }

    return( invisible(TRUE) )
    }





##############################      deadwood below        ######################################


#   arguments:
#
#   N       dimension of the cube, must be a positive integer
#   crange  range for the count of +1/2s for the edges, does not affect the vertex
#
#   returns list with components:
#       N       the input N
#       vertex  (N*(N-1)+2)x2 integer matrix with code for the vertex
#               the 1st int is the # of ones, and the 2nd is the starting position
#       edge    (2N*(N-2) + 2N) x 2 integer matrix with starting and ending index (in vertex) of the edges
#               this number applies only when crange=c(0L,N)
#

trans2subcomplex_old <- function( N, crange=c(0L,N) )
    {
    N   = as.integer(N)

    ok  = length(N)==1  &&  0<N
    if( ! ok )
        {
        log_level( ERROR, "N is invalid." )
        return(NULL)
        }

    ok  = is.numeric(crange)  &&  length(crange)==2  &&  0<=crange[1]  && crange[1]<crange[2]  && crange[2]<=N
    if( ! ok )
        {
        log_level( ERROR, "crange is invalid." )
        return(NULL)
        }

    vertex  = matrix( 0L, nrow=N*(N-1)+2, ncol=2 )

    colnames(vertex)    = c( "count", "start" )

    vertex[ 1, ]    = c( 0L, NA_integer_ )   #   south pole

    k   = 2L

    for( i in seq_len(N-1) )
        {
        vertex[ k:(k+N-1L), ]   = cbind( i, 1L:N )
        k   = k + N
        }

    vertex[ nrow(vertex), ] = c( N, NA_integer_ )   #   north pole




    out = list()


    out$N       = N
    out$vertex  = vertex

    if( N == 1 )
        {
        #   trivial case, only 1 edge, and only 2 vertices,  the 1-cube
        out$edge    = matrix( 1:2, nrow=1 )
        return(out)
        }

    midrange    = pmin( pmax(crange,1L), N-1L )

    edges   = 0L
    if( crange[1] == 0 )    edges   = edges + N

    seque   = midrange[1] + seq_len( max( diff(midrange), 0 )  ) - 1L    #; print(seque)

    edges   = edges + 2*N*length(seque)

    if( crange[2] == N )    edges   = edges + N

    edge    = matrix( 0L, nrow=edges, ncol=2 )

    i   = 1L
    if( crange[1] == 0 )
        {
        #   add the "cap" at the south pole
        for( start in 1:N )
            {
            edge[i, ]   = c( vtxidxfromcode(N,0L,NA), vtxidxfromcode(N,1L,start) )
            i = i+1L
            }
        }

    start_prev  = c(N,1L:(N-1L))      # lookup table

    for( count in seque )
        {
        for( start in 1:N )
            {
            #   find index of current vertex
            idx = vtxidxfromcode(N,count,start)

            #   add 1 on the left
            sprev   = start_prev[start]

            edge[i, ]   = c( idx, vtxidxfromcode(N,count+1L,sprev) )
            i   = i+1L

            #   add 1 on the right
            edge[i, ]   = c( idx, vtxidxfromcode(N,count+1L,start) )
            i   = i+1L
            }
        }

    if( crange[2] == N )
        {
        #   add the "cap" at the north pole
        for( start in 1:N )
            {
            edge[ i, ]    = c( vtxidxfromcode(N,N-1L,start), vtxidxfromcode(N,N,NA) )
            i   = i+1L
            }
        }

    out$edge    = edge

    return(out)
    }

#   parameters:
#   N       the dimension of the cube, must be a positive integer
#   count   the number of 1s in the vertex
#   start   the starting index of the run of 1s, 1-based
#
#   returns the row index in vertex[], as returned by trans2subcomplex()

vtxidxfromcode <- function( N, count, start )
    {
    if( count == 0 )    return( 1L )

    if( count == N )    return( N*(N-1L) + 2L )

    return( N*(count-1L) + start + 1L )
    }




#   parameters:
#   N       the dimension of the cube, must be a positive integer
#   count   integer M-vector with the number of 1s in the vertex
#   start   integer M-vector with the starting index of the run of 1s, 1-based
#
#   returns an MxN matrix where the i'th row is the desired vertex of the N-cube, [-1/2,1/2]^N

vertexfromcode_old <- function( N, count, start )
    {
    M   = length(count)
    if( length(start) != M )    return(NULL)

    #   make wrap-around lookup table
    wrap2   = c( 1:N, 1:N )

    out = matrix( -1/2, nrow=M, ncol=N )

    for( i in 1:M )
        {
        if( count[i] == 0 ) next   # south pole, all -1/2  already

        if( count[i] == N ) {
            #   north pole, all 1/2
            out[ i, ]   = 1/2
            next
        }

        jvec    = wrap2[ start[i]:(start[i] + count[i]-1L) ]
        out[ i, jvec ]  = 1/2
        }

    return(out)
    }




#   matgen      M x N matrix of N generators, defining a 2-transition subcomplex in R^M
#   centered    center the output coordinates
#
#   returns data.frame with N*(N-1)/2 rows and these columns
#       idxpair     integer matrix with 2 columns i and j with 1 <= i < j <= n
#       center      real matrix with 3 columns containing the corresponding facet center
#
#   the row order is the same as the column order returned by allcrossproducts()

allfacetcenters2trans_oldold <- function( matgen, centered=TRUE )
    {
    ok  = is.numeric(matgen)  &&  is.matrix(matgen)
    if( ! ok )  return(NULL)

    n   = ncol(matgen)
    m   = nrow(matgen)

    gensum  = .Call( C_cumsumMatrix, matgen, 2L )

    idxpair = .Call( C_allpairs, n )
    colnames(idxpair)   = c('i','j')

    facets  = nrow(idxpair)  #= ( n*(n-1L) ) / 2

    #   center is loaded with the *uncentered* coords
    center  = matrix( 0, facets, m )

    for( k in 1:facets )
        {
        i   = idxpair[k,1]
        j   = idxpair[k,2]

        v = ( matgen[ ,i] + matgen[ ,j] ) / 2  #+  (gensum[ , j-1L ] - gensum[ , i])

        if( i < j-1L )
            v = v + (gensum[ , j-1L ] - gensum[ , i])

        center[k, ] = v
        }

    if( centered )
        {
        centerpoly  = gensum[ ,n] / 2
        center      = .Call( C_sumMatVec, center, -centerpoly, 1L ) # translate to the *centered* polyhedron
        }

    out         = data.frame( row.names=1:facets )
    out$idxpair = idxpair
    out$center  = center

    return( out )
    }


allfacetcenters2trans_old  <- function( matgen, centered=TRUE )
    {
    ok  = is.numeric(matgen)  &&  is.matrix(matgen)
    if( ! ok )  return(NULL)

    n   = ncol(matgen)
    m   = nrow(matgen)

    facets    = ( n*(n-1L) ) / 2

    gensum  = .Call( C_cumsumMatrix, matgen, 2L )

    idxpair = matrix( 0L, facets, 2 )
    colnames(idxpair)   = c('i','j')

    #   center is loaded with the *uncentered* coords
    center  = matrix( 0, facets, m )

    k   = 1L

    for( i in 1:(n-1) )
        {
        for( j in (i+1):n )
            {
            idxpair[k, ]    = c(i,j)

            v   = ( matgen[ ,i] + matgen[ ,j] ) / 2

            if( i < j-1L )
                v = v + gensum[ , j-1L ] - gensum[ , i ]

            center[k, ] = v

            k   = k + 1L
            }
        }

    if( centered )
        {
        centerzono  = gensum[ ,n] / 2
        center      = .Call( C_sumMatVec, center, -centerzono, 1L )     # translate to the *centered* zonohedron
        }

    out         = data.frame( row.names=1:facets )
    out$idxpair = idxpair
    out$center  = center

    return( out )
    }

