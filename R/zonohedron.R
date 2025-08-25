#
#   zonohedron is a 3-dimensional zonotope
#
#   implemented as a list with items:
#       matroid     the matroid for it, which includes the generating matrix, etc.
#       center      a 3-vector; the center of the zonohedron
#       facet       a data.frame with a row for each facet-pair of the zonohedron,
#                   but data for only one of the facets is stored in this data.frame,
#                   and the antipodal facet is derived from this one.
#               The length and order of facet[] is the same as the order of the hyperplane list
#               in the (simplified) matroid.
#               The data.frame has these columns:
#                   center  the center of the facet, in the centered zonohedron, the other facet is antipodal (opposite sign).
#                           center is computed one belt at a time; for the first facet in the belt it is computed by maximizing wrt unit normal.
#                           the remaining facets are computed by traversing the belt, and updating the center as we go along
#                   normal  outward-pointing unit normal, the antipodal facet has the opposite normal
#                   beta    equation of the slab is  -beta <= <x,normal> <= beta.  We always have beta>0.
#                   sign    +1 or -1.  It's the difference between the facet normal and the crossproduct coming from the matroid hyperplane.
#       facet0      indexes (from facet) of those facets that contain 0
#       beltlist    the i'th entry in this list is 1/2 of the belt of the i'th generator, as a vector of facet indexes
#       zonogon     list of zonogons for the non-trivial facets. Indexing is the same as simplify(x$matroid)$hyperplane
#       frame3x2    list of 3x2 matrices for non-trivial facets.  frame3x2[[i]] maps the facet plane to the zonogon[[i]] in the XY-plane.
#       signtile    list of integer vectors, all +1 or -1, for the non-trivial facets.
#                   It's the difference between the crossprod for the pgram and the facet normal
#       zonoseg     list of zonosegs for the multiple groups.  Only present if there are multiple groups


#   zonohedron constructor
#
#   mat     a numeric matrix with 3 rows
#   e0      threshold for a column vector to be considered 0, used when nrow(x) >= 1
#   e1      threshold for codirectionality,     used when nrow(x) >= 2
#   e2      threshold for hyperplane normals,   used when nrow(x) == 3
#   ground  ground set, an integer vector in increasing order and length(ground) = ncol(x)


zonohedron <- function( mat, e0=0, e1=1.e-6, e2=1.e-10, ground=NULL )
    {
    timermain = createtimer()

    ok  = is.matrix(mat)  &&  is.numeric(mat)  &&  nrow(mat)==3  &&  3<=ncol(mat)
    if( ! ok )
        {
        log_level( ERROR,  "mat is invalid." )
        return(NULL)
        }

    matroid = matroid( mat, e0=e0, e1=e1, e2=e2, ground=ground )
    if( is.null(matroid) )
        return(NULL)

    timermain   = updatetimer(timermain)
    timematroid = timermain$elapsed

    out = list()

    class(out)  = c( "zonohedron", "zonotope", class(out) )

    #   copy non-special attributes from matrix to matroid
    for( a in names(attributes(mat)) )
        {
        if( ! ( a %in% c('class','dim','dimnames','names') ) )     attr(out,a) = attr(mat,a)
        }

    out$matroid = matroid

    out$center  = 0.5 * rowSums( getmatrix(matroid) )

    matsimple   = getsimplified(matroid)


    #   get the simplified generators, and their ground set
    matgen  = getmatrix( matsimple )
    gndgen  = getground( matsimple )

    numgen  = length(gndgen)  # number of generators


    #   make lookup table from ground to column index
    idxfromground   = integer( gndgen[numgen] )
    idxfromground[ gndgen ] = 1:numgen


    hyperplane  = gethyperplane( matsimple )
    numhypers   = length(hyperplane)

    #   2 variables for the facet centers
    center  = matrix( NA_real_, numhypers, 3 )
    delta   = double( numhypers )

    beltmat    = .Call( C_beltmatrix, hyperplane, gndgen )
    if( is.null(beltmat) ) return(NULL)


    #   facet0 is the vector of hyperplane/facet indexes that meet 0, for the simplified generators
    facet0  = integer(0)


    hyperplanes     = rep( NA_integer_, numgen )

    timebuild       = rep( NA_real_, numgen )
    timesort        = rep( NA_real_, numgen )
    timecenter1     = rep( NA_real_, numgen )
    timebeltdata    = rep( NA_real_, numgen )
    timecentersassign   = rep( NA_real_, numgen )
    timetotal       = rep( NA_real_, numgen )


    timermain   = updatetimer(timermain)
    timeprep    = timermain$elapsed

    #   edge0 is the vector of simplified generators with an edge that meets 0
    #   start with a logical mask and convert to vector later, using which()
    #edge0   = logical(numgen)

    #   beltlist[[i]] holds the 'zone' of the i'th generator, as a vector of hyperplane indexes
    beltlist    = vector( numgen, mode='list' )

    for( i in 1:numgen )
        {
        timerbelt   = createtimer()

        gen         = matgen[ ,i]               # gen is the i'th generator

        #   get the indexes of all hyperplanes in the i'th belt
        hyperidx    = beltmat[i, ]
        hyperidx    = hyperidx[ 0<hyperidx ]    # trim trailing 0s
        #cat( "------------------  gen=", i, "   hyperidx=", hyperidx, '\n' )

        m   = length(hyperidx)

        hyperplanes[i]  = m


        #   get normals to the facet-pairs that have i'th generator as an edge
        #   all these normals are orthogonal to gen, and so lie on a great circle
        normalfacet = getnormal( matsimple, hyperidx ) #; print(normalfacet)

        normalfacet = t(normalfacet)

        #   make a scratch data.frame suitable for ordering these points on the great circle
        df  = data.frame( hyperidx=rep(hyperidx,2), sign=c(rep(1L,m),rep(-1L,m)) )

        #   make the facet normal for the complete belt, including both facets in a pair
        #   but not in the correct order
        df$normal  = rbind(normalfacet,-normalfacet)
        #print( str(df) )

        if( FALSE  &&  ! all( is.finite(df$normal) ) )
            {
            log_level( FATAL, "Internal Error. in df$normal." ) ;   #  print( df$normal )
            return(NULL)
            }

        timerbelt       = updatetimer(timerbelt)
        timebuild[i]    = timerbelt$elapsed



        #   unitize the generator
        unit    = gen / sqrt( sum(gen^2) )

        #   put the true facet normals in order around the great circle
        perm    = orderoncircle( df$normal, unit )     #; print(perm) ; print( df[perm, ] )

        #   reorder and only keep the first m normals, every hyperplane will be represented exactly once
        df  = df[ perm[1:m], ]

        #   record the sorted hyperplane indexes for later use
        beltlist[[i]]   = df$hyperidx


        #   find the first hyperplane with positive sign
        #ifirst  = which( 0 < df$sign )[1]

        # keep only m of these hyperplanes, starting with ifirst
        #df  = df[ ifirst:(ifirst+m-1L), ]   #; print(df)

        timerbelt   = updatetimer(timerbelt)
        timesort[i] = timerbelt$elapsed


        #   find the center of hyperplane df$hyperidx[1]
        hyper           = matsimple$hyperplane[[ df$hyperidx[1] ]]
        generatoridx    = match( hyper, gndgen )

        pcube   = as.double(df$normal[1, ] %*% matgen)

        # force the vertex at generators to be exactly 0
        pcube[ generatoridx ] = 0
        pcube  = 0.5 * sign( pcube )

        # center1 = matgen %*% pcube        don't remember what center1 was used for

        timerbelt       = updatetimer(timerbelt)
        timecenter1[i]  = timerbelt$elapsed


        res = getbeltdata( matsimple, df$hyperidx, pcube, gndgen[i], df$normal )
        if( is.null(res) )  return(NULL)

        radmat      = res[[1]]

        centermat   = res[[2]]

        #print( res[[3]] )  # res[[2]]  is a logical mask
        if( any( res[[3]] ) )
            {
            facet0      = c( facet0, df$hyperidx[ res[[3]] ] )
            #edge0[i]    = TRUE
            }

        if( FALSE )
            {
            #   verify that radmat has the right direction, using facet normal
            for( k in 1:m )
                {
                #   get the k'th hyperplane in this belt
                # hyper       = hyperplane[[ df$hyperidx[k] ]]

                #   the "diameter vector" of this facet is a signed linear combination of the "other" generators
                #D[ ,k]  = df$sign[k] * getdiameter( matsimple, df$hyperidx[k], gndgen[i] )

                #   do a sign test
                d2  = sum( crossproduct(radmat[ ,k], gen) * df$normal[k, ] )

                if( d2 <= 0 )
                    {
                    log_level( FATAL, "Internal sign error. k=%d.  d2=%g <= 0.", k, d2 )
                    return(NULL)
                    }
                }
            }

        timerbelt       = updatetimer(timerbelt)
        timebeltdata[i] = timerbelt$elapsed


        #   load centermat into center.
        #   since each facet is contained in 2 or more belts,
        #   centers will be assigned 2 or more times.
        #   The array delta[] records the tiny differences.  It is modified in place.
        #   the vector df$sign is replicated to all 3 columns
        #center[ df$hyperidx, ] =  df$sign * t(centermat)
        res = .Call( C_multicopy, center, delta, df$sign * t(centermat), df$hyperidx )
        if( is.null(res) )
            {
            log_level( FATAL, "Internal Error.  C_multicopy() fail.  i=%d.", i )
            return(NULL)
            }

        timerbelt               = updatetimer(timerbelt)
        timecentersassign[i]    = timerbelt$elapsed
        timetotal[i]            = timerbelt$total
        }

    timermain   = updatetimer(timermain)
    timebelts   = timermain$elapsed


    #   test that all facet centers were assigned
    if( any( is.na(center) ) )
        {
        log_level( ERROR, "some facet centers not assigned." )
        return(NULL)
        }

    #   test that facet centers have negligible disagreement
    centernorm  = .rowSums( abs(center), nrow(center), ncol(center) )
    mask    = (centernorm == 0)
    if( any(mask) )
        {
        log_level( ERROR, "%d facet centers are 0.", sum(mask) )
        return(NULL)
        }
    deltarel    = delta / centernorm    #; cat( "deltarel=", deltarel, '\n' )
    tol = 5.e-12
    mask    = tol < deltarel
    if( any(mask) )
        {
        log_level( WARN, "%d facet centers have disagreement > %g.", sum(mask), tol )
        }

    normal  = t( getnormal( matsimple, NULL ) )     # get ALL the normal vectors

    if( TRUE  &&  any( out$center!=0 ) )
        {
        #   modify center and normal so the *chosen* facet of the pair is closer to 0 (uncentered)
        #   beta[] remains unchanged
        #   signvec = -as.numeric( sign( center %*% out$center ) )
        signvec = -sign( normal %*% out$center )

        signvec[ signvec==0 ]   = 1    # change any 0s to 1s,  should be extremely rare

        #center  = signvec * center      # uses recycling rule so signvec is multiplied by all columns
        #normal  = signvec * normal      # uses recycling rule so signvec is multiplied by all columns

        .Call( C_timesEqual, center, signvec, 2L )     # multiply in place
        .Call( C_timesEqual, normal, signvec, 2L )     # multiply in place
        }
    else
        {
        #   do not modify center and normal
        signvec = rep( 1L, nrow(center) )
        }


    #   build the facet data.frame and add to output list
    facet           = data.frame( row.names=1:numhypers )
    facet$center    = center
    facet$normal    = normal
    facet$beta      = .rowSums( center*normal, nrow(center), ncol(center) )
    facet$sign      = as.integer( signvec )

    out$facet   = facet


    #   build the facet0 data.frame and (possibly) modify it for the original generators
    facet0  = sort( unique(facet0) )   # unique(facet0)
    #cat( "simplified facet0=", facet0, '\n' )

    if( FALSE )
        {
        #   convert from logical vector to indexes
        edge0   = which( edge0 )

        if( length(edge0) != length(facet0) )
            {
            log_level( FATAL, "Internal Error.  length(edge0) = %d  !=  %d = length(facet0).",
                            length(edge0), length(facet0) )
            return(NULL)
            }
        }


    #   facet0 is now for the simplified generators.
    #   modify facet0[] for the original generators, by removing some of them.
    #   get the "mixed direction" generators, which are the generator minors that are non-zero
    colidx  = getmixed( matroid )   # column indexes of the simplified matrix, but pass the original matroid

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
        }

    out$facet0  = facet0

    #out$edge0   = edge0

    out$beltlist    = beltlist

    timermain   = updatetimer(timermain)
    timefacets  = timermain$elapsed

    #   for each non-trivial facet  (non-parallelogram) compute a zonogon
    lenvec      = lengths( hyperplane )
    idxnontriv  = which( 2 < lenvec )
    hypersnt    = length(idxnontriv)
    if( 0 < hypersnt )
        {
        if( ! all( idxnontriv == 1:hypersnt ) )
            {
            log_level( FATAL, "Internal Error.  %d nontrivial hyperplanes are not contiguous.", hypersnt )
            return(NULL)
            }

        frame3x2    = vector( hypersnt, mode='list' )
        zonogon     = vector( hypersnt, mode='list' )
        signtile    = vector( hypersnt, mode='list' )

        for( k in 1:hypersnt )
            {
            #   get the facet normal and make frame from it
            #normal  = facet$normal[ k, ]

            frame3x2[[k]]   = frame3x2fun( facet$normal[ k, ] )

            #   get the column indexes for the hyperplane
            idxhyper    = idxfromground[ hyperplane[[k]] ]

            #   extract only the generators for this face
            matsub  = matgen[ , idxhyper ]

            #   rotate generators to dimension 2
            matsub2 = t( frame3x2[[k]] )  %*%  matsub

            #   construct the non-trivial zonogon
            zonogon[[k]]    = zonogon( matsub2, e0=0, e1=0, ground=gndgen[idxhyper] )
            if( is.null(zonogon[[k]]) )
                {
                log_level( FATAL, "Internal Error.  Cannot construct zonogon %d.", k )
                return(NULL)
                }

            #   each element in signtile[[k]] is +1 or -1, depending on whether the crossproduct
            #   agrees with the chosen facet normal, or is its negative
            #   get all crossprods for this non-trivial facet
            #   NOTE:  if the facet is trivial, the chosen normal *is* the cross-product
            idx = translateallpairs( hyperplane[[k]], getground(matsimple) )

            crossprods      = matsimple$crossprods[ , idx ]

            signtile[[k]]   = as.integer( sign( facet$normal[ k, ] %*% crossprods ) )
            }

        out$zonogon     = zonogon
        out$frame3x2    = frame3x2
        out$signtile    = signtile
        }

    #   for each multiple group, compute the corresponding zonoseg
    nummultiples    = length(matroid$multiple)

    if( 0 < nummultiples )
        {
        out$zonoseg     = makezonoseglist( matroid )
        }



    timermain       = updatetimer(timermain)
    timezonogons    = timermain$elapsed

    timeall     = timermain$total

    if( FALSE  )
        {
        perf    = data.frame( point=gndgen, hyperplanes=hyperplanes )
        perf$build          = timebuild * 1000
        perf$sort           = timesort * 1000
        perf$center1        = timecenter1 * 1000
        perf$beltdata       = timebeltdata * 1000
        perf$beltdatamean   = perf$beltdata/hyperplanes
        perf$centersassign  = timecentersassign * 1000
        perf$total          = timetotal * 1000
        print( perf )
        }

    if( FALSE )
        {
        cat( "matroid:      ", timematroid * 1000, " msec\n" )
        cat( "preparation:  ", timeprep * 1000, " msec\n" )
        cat( "belts:        ", timebelts * 1000, " msec.  ", "belts=", numgen, "\n" )
        cat( "facets:       ", timefacets * 1000, " msec\n" )
        cat( "zonogons:     ", timezonogons * 1000, " msec\n" )
        cat( "total:        ", timeall * 1000, " msec\n" )
        }

    return( out )
    }

#   n       a positive integer, so the step size on the circle is 2*pi/n
#   m       number of points to compute, starting at 1
#   height  height (of "white point") when m==n

polarzonohedron <- function( n, m=n, height=pi, ground=NULL )
    {
    if( is.null(ground) )
        ground  =  1:m
    else if( length(ground) != m )
        {
        log_level( ERROR, "ground is invalid, because the length is incorrect." )
        return(NULL)
        }

    if( m<3 || n<m )
        {
        log_level( ERROR, "m=%d is invalid.", m )
        return(NULL)
        }

    u   = (0:(m-1)) / n

    mat = t( tocircle(2*pi*u) )
    mat = rbind( mat, 1 ) * height/n

    return( zonohedron(mat,e0=0,e1=0,e2=0,ground=ground) )
    }





#   mat     2xN matrix giving N points in the plane, N>=3  Assumed to be already simplified, this is verified.
#           the order of the column generators does not affect the returned zonohedron

#   returns zonohedron with matroid simple and uniform, these are verified.
#           As customary, only half of the facets are returned, and the antipodal ones are comitted.
#           Only the "lower" facets are returned, these project down to a regular tiling of the zonogon defined by mat.
#           These "lower" facets are the ones visible from *below*, on the xy-plane.
#           out$facet is modified so all the facets have normals with negative Z-component

liftedzonohedron <- function( mat, e0=0, e1=0, e2=0, ground=NULL )
    {
    ok  = is.matrix(mat)  &&  is.numeric(mat)  &&  nrow(mat)==2  &&  3<=ncol(mat)
    if( ! ok )
        {
        log_level( ERROR, "argument mat is invalid." )
        return(NULL)
        }

    #   'lift' generators to a cone with vertex at 0
    mat3    = rbind( mat, sqrt(colSums(mat^2)) )

    out     = zonohedron( mat3, e0=e0, e1=e1, e2=e2, ground=ground )

    if( is.null(out) )  return(NULL)

    if( ! is_simple(out$matroid) )
        {
        log_level( ERROR, "matroid is not simple." )
        return(NULL)
        }

    if( ! is_uniform(out$matroid) )
        {
        log_level( ERROR, "matroid is not uniform." )
        return(NULL)
        }

    #   force all facet normals to point "down", so all facets are on the "lower" side
    signz   = sign( out$facet$normal[ ,3] )

    if( any(signz==0) )
        {
        log_level( ERROR, "%d facets have horizontal facet normal.", sum(signz==0) )
        return(NULL)
        }

    maskpos = (0 < signz)
    if( any(maskpos) )
        {
        #   reverse normal and center, so all facets in the data.frame are on the "lower" side
        out$facet$normal[maskpos, ] = -out$facet$normal[maskpos, ]
        out$facet$center[maskpos, ] = -out$facet$center[maskpos, ]
        }

    return( out )
    }



#   n       a positive integer, so the step size on the circle is 2*pi/n
#   m       number of points to compute, starting at 1
#   axis    extrusion axis, z must be nonzero

regularprism <- function( n, m=n, axis=c(0,0,1), ground=NULL )
    {
    if( is.null(ground) )
        ground  = 1:(m+1)
    else if( length(ground) != m+1 )
        {
        log_level( ERROR, "ground is invalid, because the length is incorrect." )
        return(NULL)
        }

    if( m<1 || n<m )
        {
        log_level( ERROR, "m=%d is invalid.", m )
        return(NULL)
        }

    ok  = is.numeric(axis) && (length(axis) %in% c(1,3))
    if( ! ok )
        {
        log_level( ERROR, "axis is invalid." )
        return(NULL)
        }

    if( length(axis) == 1 )
        axis    = c(0,0,axis)

    if( axis[3] == 0 )
        {
        log_level( ERROR, "axis is invalid, because its z coordinate is 0." )
        return(NULL)
        }

    u   = (0:(m-1)) / n

    mat = t( tocircle(2*pi*u) )     # m columns
    mat = rbind( mat, 0 )           # m columns
    mat = cbind( mat, axis )        # m+1 columns

    return( zonohedron(mat,ground=ground) )
    }


quasicube <- function( count=c(1,1,1) )
    {
    count   = as.integer(count)

    if( length(count) == 1 )
        count   = rep( count, 3 )
    else if( length(count) == 2 )
        count   = c( count, 0L )

    ok  = length(count)==3  && all(0<=count)  &&  sum( (count==0) <=1 )
    if( ! ok )
        {
        log_level( ERROR, "vector count=%d,%d,%d is invalid", count[1], count[2], count[3] )
        return(NULL)
        }

    n   = count[1]
    if( 0 < n )
        {
        #   xy plane
        u   = seq( 0, by=2*pi/n, length.out=n )
        mat1    = t( cbind( tocircle(u), 0 ) )
        }
    else
        mat1    = matrix( 0, nrow=3, ncol=0 )


    n   = count[2]
    if( 0 < n )
        {
        #   yz plane
        u   = seq( 0, by=2*pi/n, length.out=n )
        mat2    = t( cbind( tocircle(u), 0 ) )
        mat2    = mat2[ c(3,1,2), , drop=F]
        }
    else
        mat2    = matrix( 0, nrow=3, ncol=0 )


    n   = count[3]
    if( 0 < n )
        {
        #   zx plane
        u   = seq( 0, by=2*pi/n, length.out=n )
        mat3    = t( cbind( tocircle(u), 0 ) )
        mat3    = mat3[ c(2,3,1), , drop=F]
        }
    else
        mat3    = matrix( 0, nrow=3, ncol=0 )

    matgen  = cbind( mat1, mat2, mat3 ) #; print(matgen)

    out = zonohedron( matgen )

    return( out )
    }


print.zonohedron  <-  function( x, trans2=FALSE, matroid=TRUE, ... )
    {
    cat( "zonohedron:\n" )

    ground  = getground( getsimplified(x$matroid) )

    fullname    = attr( x, "fullname" )
    if( ! is.null(fullname) )
        cat( "fullname:                         ", fullname, "\n" )

    gens    = length( getground(x$matroid) )
    cat( "generators (original):            ", gens, '\n' )

    gens   = length( x$matroid$multiple )
    cat( "generators with multiples:        ", gens, '\n' )

    if( 0 < gens )
        {
        colidx  = getmixed(x$matroid)   #;which( 0 < rowSums( abs(x$matroid$multiplesupp$minor) ) )
        gens    = length(colidx)
        cat( "generators with mixed-directions: ", gens,  '   {', ground[colidx], '}\n' )
        }

    gens    = length( getground(getsimplified(x$matroid)) )
    cat( "generators (simplified):          ", gens, '\n' )

    facets  = 2*nrow(x$facet)
    cat( "number of facets:                  ", facets, '  [', facets/2, " antipodal facet-pairs]\n", sep='' )
    cat( "facets that contain 0:            ", length(x$facet0), "   {", x$facet0, "}\n" )

    res = getmetrics( x )
    cat( "number of edges:                  ", res$edges, "\n"  )
    # cat( "edges that contain 0:             ", length(x$edge0), "   {", ground[x$edge0], "}\n" )


    cat( "center:                           ", x$center, '\n' )
    cat( "pointed:                          ", is_pointed(x), '\n' )
    cat( "salient:                          ", is_salient(x), '\n' )

    cat( "area:                             ", res$area, '\n' )
    cat( "volume:                           ", res$volume, '\n' )

    flush.console()

    if( trans2 )
        {
        extra   = getmetrics2trans( x )

        cat( '\n' )
        cat( "2-Transition Polyhedron Metrics:\n" )

        cat( "orientation:                      ", extra$signsurf, '\n' )

        mess    = sprintf( " [signcounts: %d + %d + %d = %d]",
            extra$signcount[1], extra$signcount[2], extra$signcount[3], extra$signcount[4] )

        cat( "strictly starshaped at center:    ", extra$starshaped, mess, '\n' )
        cat( "injective:                        ", extra$injective, '\n' )

        if(  is.finite(extra$starshaped)  &&  extra$starshaped )
            {
            deficient   = extra$pgramdf$deficient

            dfacets = 2 * sum( deficient )
            cat( "deficient parallelograms:         ", dfacets, "   [fraction=", dfacets/facets, ']\n' )
            cat( "deficient area:                   ", extra$areadeficient, "   [fraction=", extra$areadeficient/res$area, ']\n' )
            cat( "deficient volume:                 ", extra$volumedeficient, "   [fraction=", extra$volumedeficient/res$volume, ']\n' )

            thickness = extra$volumedeficient / extra$areadeficient
            cat( "deficient thickness (mean):       ", thickness, '\n' )

            #   compute transitions
            pgramdef    = extra$pgramdf[ deficient, ]

            dat         = boundarypgramdata( x, pgramdef$gndpair  )
            #dat$gndpair = pgramdef$gndpair
            dat$area    = pgramdef$area

            #   breakdown the deficient parallelograms by transitions
            tunique = sort( unique( dat$transitions )  )    #; print( tunique )
            m   = length(tunique)

            pgcount = integer(m)
            area    = numeric(m)
            example = character(m)

            for( k in 1:m )
                {
                datsub      = dat[ dat$transitions == tunique[k], ]
                pgcount[k]  = nrow( datsub )
                area[k]     = sum( datsub$area )

                i  = which.max( datsub$area )
                gndpair     = datsub$gndpair[i, ]
                example[k]  = sprintf( "{%d,%d}", gndpair[1], gndpair[2] )
                }

            dftemp  = data.frame( transitions=c(tunique, "Totals") )
            dftemp$parallelograms   = c( 2L*pgcount, 2L*sum(pgcount) )
            dftemp$area             = c( 2*area, 2*sum(area) )
            dftemp$example          = c( example, '' )
            cat( '\n' )
            print( dftemp, row.names=FALSE )

            cat( '\n' )
            edges_total     = 2L * nrow(extra$anglesDH)
            edges_concave   = 2L * sum( extra$anglesDH$angle < 0 )
            convex  = edges_concave==0

            cat( "total edges:              ", edges_total, '\n' )
            cat( "concave edges:            ", edges_concave, "   [fraction=", edges_concave/edges_total, ']\n' )

            if( ! convex )
                {
                theta       = pi + extra$anglesDH$angle
                imin        = which.min( theta )    #; print( imin )

                thetamin    = theta[imin]

                dfsub       = extra$anglesDH[imin, ]

                pivot   = ground[ dfsub$pivot[1] ]
                wing1   = ground[ dfsub$wing1[1,1] ]
                wing2   = ground[ dfsub$wing2[1,1] ]

                mess    = sprintf( "at edge common to parallelograms {%g,%g} and {%g,%g}",
                                            pivot, wing1, pivot, wing2 )

                cat( "minimum external angle:   ", thetamin, 'radians, ', mess, '\n' )
                }
            cat( "convex surface:           ", convex, '\n' )
            }
        else
            {
            cat( '\n' )
            cat( "****  Cannot print more 2-Transition metrics, since the surface is not strictly starshaped at the center.  ****\n" )
            }
        }

    if( matroid )
        {
        cat( '\n' )
        cat( "matroid:\n" )
        print( x$matroid )
        }

    return( invisible(TRUE) )
    }



#   returns a list with items:
#       vertices    number of vertices, computed from edges and facets using Euler characteristic
#       edges       number of edges
#       facets      number of facets
#       area        surface area
#       volume      volume


getmetrics.zonohedron <- function( x )
    {
    matsimple   = getsimplified( x$matroid )

    facets      = 2L * nrow(x$facet)

    edges       = 2L * sum( lengths(gethyperplane(matsimple)) )

    vertices    = edges - facets + 2L

    #   matrix crossprods[,] was already computed when the matroid was, but was then unitized
    #   NOTE:  this matrix crossprodsraw[,] is raw and nonunitized
    crossprodsraw   = allcrossproducts( getmatrix(matsimple) )

    areavec = sqrt( .colSums( crossprodsraw^2, nrow(crossprodsraw), ncol(crossprodsraw) ) )

    out = list()

    out$vertices    = vertices
    out$edges       = edges
    out$facets      = facets

    out$area    = 2*sum(areavec)

    #   matsimple$hyperplaneidx is a LUT from crossprod index to hyperplane index
    #   so for non-trivial hyperplanes, the beta is duplicated many times
    betavec     = x$facet$beta[ matsimple$hyperplaneidx ]   #; print(betavec)

    if( length(areavec) != length(betavec) )
        {
        log_level( FATAL, "Internal Error.   length(areavec)=%d != %d=length(betavec).",
                 length(areavec), length(betavec) )
        #print( matsimple$hyperplaneidx )
        #print( str(betavec) )
        return(out)
        }

    out$volume  = (2/3) * sum( areavec * betavec )

#    if( trans2 )
#        {
#        extra   = getmetrics2trans( x, tol )
#        if( is.null(extra) )  return(NULL)
#
#        out$trans2  = extra
#        }

    return(out)
    }




#   x   a zonohedron
#   i   index of a hyperplane of the simplified matroid of x

gettilecenters3D  <- function( x, i )
    {
    if( i <= length(x$zonogon) )
        {
        #   non-trivial hyperplane, many tiles
        zono    = x$zonogon[[i]]

        #   find center of all pgram tiles in 3D
        center3D    = zono$tilingdata$center %*% t( x$frame3x2[[i]] )

        #   translate each pgram center to the zonogon facet in 3D
        center3D    = .Call( C_sumMatVec, center3D, x$facet$center[i, ], 1L )
        }
    else
        {
        #   non-trivial hyperplane, 2 generators, only 1 tile = the pgram itself
        center3D    = x$facet$center[i, , drop=FALSE]
        }

    return( center3D )
    }





if( FALSE )
{

#   x   a zonohedron object
#   p   an Mx3 matrix, etc.
#
#   value   see inside_zonotope()

inside.zonohedron <- function( x, p )
    {
    p   = prepareNxM( p, 3 )
    if( is.null(p) )    return(NULL)

    #   translate p to the centered zonohedron
    gcentered   = p - matrix( x$center, nrow(p), 3, byrow=TRUE ) #; print(gcentered)

    hg  = tcrossprod( x$facet$normal, gcentered )    #; print( str(hg) )

    out = inside_zonotope( x, p, hg )

    return( out )
    }
}



#   x       a zonohedron object
#   ...     possibly even *more* zonohedron objects !!
#
#   return  a data.frame with a row for each zonohedron, with important metrics.  N rows

summary.zonohedron <- function( object, ... )
    {
    #   combine all the zonohedron objects into a single list
    zlist = c( list(object), list(...) ) #; print( zlist )

    summary_from_zlist( zlist )
    }


#   zlist   a list of zonohedron object
#   full    if TRUE, then include area and volume and pointed columns
#
#   return  a data.frame with a row for each zonohedron, with important metrics.  N rows

summary_from_zlist <- function( zlist, full=TRUE )
    {
    ok  = sapply( zlist, inherits, what="zonohedron" )      #; print(ok)

    if( ! all(ok) )
        {
        log_level( ERROR, "%d of the %d objects in the list are not zonohedra.", sum(!ok), length(ok) )
        return(NULL)
        }

    n   = length(zlist)

    shortname   = rep( NA_character_, n )
    fullname    = rep( NA_character_, n )
    generators  = rep( NA_integer_, n )
    #center      = matrix( NA_real_, nrow=n, ncol=3 )
    vertices    = rep( NA_integer_, n )
    edges       = rep( NA_integer_, n )
    facets      = rep( NA_integer_, n )
    area        = rep( NA_real_, n )
    volume      = rep( NA_real_, n )
    pointed     = rep( NA, n )

    for( k in 1:n )
        {
        zono    = zlist[[k]]

        sn  = attr(zono,"shortname")
        if( ! is.null(sn) ) shortname[k] = sn

        fn  = attr(zono,"fullname")
        if( ! is.null(fn) ) fullname[k] = fn

        generators[k]   = ncol( getmatrix( zono$matroid ) )

        #   center[k, ] = zono$center

        mets        = getmetrics(zono)

        vertices[k] = mets$vertices
        edges[k]    = mets$edges
        facets[k]   = mets$facets

        area[k]     = mets$area
        volume[k]   = mets$volume

        pointed[k]  = is_pointed( zono )
        }

    rnames  = names(zlist)
    if( is.null(rnames)  ||  anyDuplicated(rnames)!=0 )
        rnames  = 1:n

    out = data.frame( row.names=rnames )

    #   out$shortname   = shortname
    out$fullname    = fullname
    out$generators  = generators
    #$out$center      = center
    out$vertices    = vertices
    out$edges       = edges
    out$facets      = facets

    if( full )
        {
        out$area    = area
        out$volume  = volume
        out$pointed = pointed
        }

    return( out )
    }


#   x           a zonohedron object
#   base        a numeric vector of length 3, the basepoint of all the rays
#               base must be in the interior of x,
#               or if x is non-negative, base can also be the black or white point on the boundary(x)
#   direction   an Nx3 matrix with directions in the rows; each direction must be non-zero
#
#   value   a dataframe with columns
#           base        given basepoint of all the rays (all the same)
#           direction   given directions of the rays
#           facetidx    idx of the facet where ray exits the zonohedron
#           sign        +1 if beta, and -1 if -beta
#           tmax        ray parameter of intersection with facet
#           point       point of intersection ot the ray and the facet
#           timetrace   time to do the trace, in sec
#

raytrace.zonohedron <- function( x, base, direction, invert=FALSE, plot=FALSE, ... )
    {
    ok  = is.numeric(base)  &&  length(base)==3  &&  all( is.finite(base) )
    if( ! ok )
        {
        log_level( ERROR, "base is invalid. It must be a numeric vector of length 3, and all entries finite." )
        return(NULL)
        }

    direction   = prepareNxM( direction, 3 )
    if( is.null(direction) )    return(NULL)

    #   translate base to the centered zonohedron
    base    = as.numeric(base)

    gcentered   = base - x$center  #; print(gcentered)

    dim(base)       = c(1,3)
    dim(gcentered)  = c(1,3)

    hg  = tcrossprod( x$facet$normal, gcentered )    #; print( str(hg) )

    # hg  = as.numeric( x$facet$normal  %*%  gcentered )      #; print( str(hg) )


    #   test whether base is black or white point, no tolerance here
    dim(gcentered)  = NULL

    blackpt = ifelse( is_salient(x),  all(gcentered == -x$center), FALSE )
    whitept = ifelse( is_salient(x),  all(gcentered ==  x$center), FALSE )

    if( blackpt || whitept )
        {
        #   get the normals for all facets that meet 0
        # normal0 is Mx3 where M is the number of these facets
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
            log_level( ERROR, "point base=(%g,%g,%g) is not in the interior of the zonohedron. distance=%g >= 0",
                                    base[1], base[2], base[3], df$distance )
            return(NULL)
            }
        }


    dim(base)       = NULL



    n   = nrow(direction)

    tmax        = rep(NA_real_,n)
    idx         = rep(NA_integer_,n)
    sign        = rep(NA_integer_,n)
    point       = matrix(NA_real_,n,3)
    timetrace   = rep(NA_real_,n)

    for( k in 1:n )
        {
        time_start  = gettime()

        v   = direction[k, ]

        if( any( is.na(v) ) )   next

        if( sum(v*v) == 0 ) next    # 0-vector

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

        if( tvec[j] <= 0 )
            {
            log_level( WARN, "Internal Warning.  direction=(%g,%g,%g) failed to intersect the boundary properly.",
                                    v[1], v[2], v[3] )
            next    # failed to intersect properly
            }

        tmax[k] = tvec[j]           # tmax[k] is positive

        idx[k]  = j  # x$facet$idx[j]

        sign[k] = as.integer( sign(hv[j]) )

        point[k, ]      = base  +  tmax[k] * v     #; print( optcentered )

        timetrace[k]    = gettime() - time_start
        }

    rnames  = rownames(direction)
    if( is.null(rnames)  ||  anyDuplicated(rnames) )   rnames = 1:n

    out = data.frame( row.names=rnames )

    out$base        = matrix( base, n, 3, byrow=TRUE )  # replicate base to all rows
    out$direction   = direction
    out$facetidx    = idx
    out$sign        = sign
    out$tmax        = tmax
    out$point       = point
    out$timetrace   = timetrace

    cnames  = colnames(base)
    if( is.null(cnames) )   cnames = colnames(direction)

    colnames(out$point)   = cnames

    if( invert )
        {
        dat = invertboundarydata( x, out )
        if( ! is.null(dat) )
            {
            out$distance    = dat$distance
            out$pcube       = dat$pcube
            out$transitions = dat$transitions
            }
        }

    if( plot )
        {
        if( ! requireNamespace( 'rgl', quietly=TRUE ) )
            log_level( WARN, "Package 'rgl' is required for plotting.  Please install it." )
        else if( rgl::cur3d() == 0 )
            log_level( WARN, "Cannot add raytrace to plot, because there is no rgl window open." )
        else
            {
            xyz = matrix( base, nrow=n, ncol=3, byrow=T )
            xyz = rbind( xyz, point )

            perm    = 1:(2*n)
            dim(perm)    = c(n,2L)
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

    return( out )
    }


#   section() compute intersection of zonohedron and plane(s)
#
#   x           a zonohedron object
#   normal      a non-zero numeric vector of length 3, the normal of all the planes
#   beta        a vector of plane-constants.  The equation of plane k is: <x,normal> = beta[k]
#
#   value   a list of length = length(beta).  Each item in the list is a list with these items:
#           beta        given plane constant of the plane
#           section     Mx3 matrix of points on the section = a polygon, in order around the boundary
#                       M=1 for a supporting plane, and if there is no intersection, then M=0 rows.

section.zonohedron <- function( x, normal, beta, tol=1.e-10, plot=FALSE, ... )
    {
    timermain = createtimer()

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

    #   unitize normal
    #normal  = normal / sqrt( sum(normal^2) )

    matsimple   = getsimplified(x$matroid)
    matgen  = getmatrix( matsimple )
    gndgen  = getground( matsimple )

    numgen  = ncol( matgen )    # so matgen is 3 x numgen

    #   for each edge, compute radius of the projection of the edge onto the line generated by normal
    normalgen   = as.numeric( normal %*% matgen )
    radiusgen   = 0.5 * abs(normalgen)              # so length(radiusgen) = numgen

    timermain   = updatetimer( timermain )
    timeedges   = timermain$elapsed


    #   for each facet/hyperplane, compute radius of the projection of the facet onto the line generated by normal
    if( FALSE )
        {
        #   make lookup table from ground to column index
        idxfromground   = integer( gndgen[numgen] )
        idxfromground[ gndgen ] = 1:numgen

        myfun   <- function( hyper )    { sum( radiusgen[ idxfromground[hyper] ] ) }

        radiusfacet = sapply( matsimple$hyperplane, myfun ) #; print(radiusfacet)
        }
    else
        {
        #   faster method in plain C
        radiusfacet = .Call( C_radiusfacet, matsimple$hyperplane, gndgen, radiusgen )
        }


    timermain   = updatetimer( timermain )
    timefacets  = timermain$elapsed

    #   numfacets is really the number of facet pairs
    numfacets   = nrow( x$facet )

    #   make matrix of *all* facet centers with 2*numfacets rows
    center      = rbind( x$facet$center, -x$facet$center )

    #   make matrix of *all* facet normals with 2*numfacets rows
    normalfacet = rbind( x$facet$normal, -x$facet$normal )



    #   make matching products of these centers and normal vector
    cn  = as.numeric( x$facet$center %*% normal )
    cn  = c( cn, -cn )

    #   compute range of normal over each facet
    cnneg   = cn - radiusfacet
    cnpos   = cn + radiusfacet

    hyperplane2 = rep( matsimple$hyperplane, 2 )    # or should it be c(matsimple$hyperplane,matsimple$hyperplane) ?

    timermain   = updatetimer( timermain )
    timecenters = timermain$elapsed


    #   translate beta to the centered zonogon
    cent.norm = sum( x$center * normal )
    betacentered   = as.numeric(beta) - cent.norm  #; print(betacentered)


    #   compute radius of the entire zonohedron, after projection onto the line generated by normal
    betamax = sum( radiusgen )
    betamin = -betamax  # by symmetry

    #cat( "centered betamin=", betamin, "   betamax=", betamax, '\n' )
    #cat( "betacentered=", betacentered, '\n' )

    res = support( x, normal )
    argmax  = res$argmax
    argmin  = 2*x$center - argmax      #-(argmax - x$center) +  cent.norm

    #cat( "betamax=", betamax, "   maxvalue - cent.norm = ", res$value - cent.norm, '\n' )


    #   find 3x2 matrix for projection to 2D
    #   this works for every section, because all are parallel
    frame3x2    = frame3x2fun( normal )

    out         = vector( length(beta), mode='list' )
    names(out)  = sprintf( "normal=%g,%g,%g. beta=%g", normal[1], normal[2], normal[3], beta )

    for( k in 1:length(beta) )
        {
        beta_k  = betacentered[k]

        if( beta_k < betamin-tol  ||  betamax+tol < beta_k )
            {
            #   plane does not intersect the zonohedron, there is no section
            df  = data.frame( row.names=character(0) )
            df$point    = matrix( 0, 0, 3 )
            out[[k]]    = df    # list( beta=beta[k],  section=matrix( 0, 0, 3 ) )
            next
            }

        if( abs(beta_k - betamax) < tol )
            {
            #   special case - only one point of intersection
            df  = data.frame( row.names=1 )
            df$point    = matrix( argmax, 1, 3 )
            out[[k]]    = df    # list( beta=beta[k],  section= matrix( argmax, 1, 3 ) )
            next
            }
        if( abs(beta_k - betamin) < tol )
            {
            #   special case - only one point of intersection
            df  = data.frame( row.names=1 )
            df$point    = matrix( argmin, 1, 3 )
            out[[k]]    = df  # = list( beta=beta[k],  section= matrix( argmin, 1, 3 ) )
            next
            }

        #   find indexes of all facets that intersect this plane
        indexvec = which( cnneg < beta_k  &  beta_k < cnpos )  #; print( indexvec )

        if( length(indexvec) == 0 )
            {
            # should not happen
            log_level( WARN, "Internal Error. length(indexvec) == 0." )
            next
            }

        #cat( "indexvec=", indexvec, '\n' )

        log_level( INFO, "beta[%d] = %g.   section has %d edges.", k, beta[k], length(indexvec) )

        #   extract the sublist of only those facets/hyperplanes that intersect the plane
        hypersub    = hyperplane2[ indexvec ]

        # in the C function, non-trivial hyperplanes are skipped and the row of section is filled with NAs
        section =  .Call( C_sectionzonohedron, hypersub, center[indexvec, ], normalfacet[indexvec, ], cn[indexvec], beta_k,
                                gndgen, normalgen, matgen, matsimple$crossprods )

        if( is.null(section) )  return(NULL)

        #   assign rownames(section)
        idxhyperplane2  = rep( 1:length(matsimple$hyperplane), 2 )

        idxhypersub    = idxhyperplane2[ indexvec ]

        antipodalmask   = (numfacets < indexvec)

        signchar    = ifelse( antipodalmask, '-', '' )

        rownames( section ) = paste( idxhypersub, signchar, sep='' )


        lenvec  = lengths(hypersub)
        if( any( 2 < lenvec ) )
            {
            #   fill in the gaps in section[]
            idxnt   = which( 2 < lenvec )

            log_level( INFO, "section k=%d has %d nontrivial facets.", k, length(idxnt) )

            for( j in 1:length(idxnt) )
                {
                idxhyper    = idxhypersub[ idxnt[j] ]
                hyper       = matsimple$hyperplane[[idxhyper]]

                #cat( "j=", j, "   idxnt[j]=", idxnt[j], "   idxhyper=", idxhyper, "  length(hyper)=", length(hyper), '\n' )

                #   for the facet center we only use the "original" non-antipodal facet
                facetcenter = x$facet$center[ idxhyper, ]

                #idx2    = indexvec[ idxnt[j] ]

                antipodal   = antipodalmask[ idxnt[j] ]

                #cat( "Antipodal = ", antipodal, '\n' )

                thesign = sign( 0.5 - antipodal  )     # so FALSE -> 1 and TRUE -> -1
                #cat( "thesign = ", thesign, '\n' )

                #print( hyper )


                #   get the zonogon and the 3x2 matrix
                zonogon = x$zonogon[[ idxhyper ]]
                A       = x$frame3x2[[ idxhyper ]]

                if( is.null(zonogon) || is.null(A) )
                    {
                    # should not happen
                    log_level( ERROR, "Internal Error. zonogon or A is NULL." )
                    return(NULL)
                    }

                #   change of coordinates from 3D to 2D
                norm2   = as.double( normal %*% A )
                temp3   = as.double( -A %*% zonogon$center )   +   facetcenter
                beta2   = thesign * beta_k - sum( temp3 * normal )
                sec2    = section( zonogon, norm2, beta2 )
                #print( sec2 )

                p2  = sec2$boundary2[1, ]

                #   now map p2 back to 3D and the centered zonohedron
                p3  = as.double( A %*% (p2 - zonogon$center) )  +  facetcenter

                section[ idxnt[j], ] = thesign * p3
                }

            #print( section )
            }


        #   check for NAs
        bad = is.na( section[ ,1] )
        if( any(bad) )
            {
            log_level( WARN, "%d points of the section (beta=%g) could not be computed, and have been set to NA.",
                        sum(bad), beta[k]  )
            }

        if( FALSE )
            {
            #   check the result
            test    = as.double( section %*% normal)    #; print(test - beta_k)
            test    =  5.e-10 < abs(test - beta_k)      #; print(test)
            if( any(test,na.rm=TRUE) )
                {
                log_level( ERROR, "Internal Error. k=%d.  %d points not on beta_k = %g plane.",
                                    k, sum(test,na.rm=TRUE), beta_k )
                print(section)
                next
                }
            }


        #   find suitable point in the interior of this section, using argmin and argmax
        s   = (beta_k - betamin) / (betamax - betamin)

        center_section = (1-s)*argmin + s*argmax  -  x$center

        #   project to 2D
        p2D = section  %*%  frame3x2

        #   subtract projection center_section, so points of section now go around the origin
        c2D = center_section  %*%  frame3x2

        res = .Call( C_plusEqual, p2D, -c2D, 1L )
        if( is.null(res) )  return(NULL)

        #   not a polygon yet, the points must be ordered by angle using atan2()
        perm    = order( atan2(p2D[ ,2],p2D[ ,1]), na.last=FALSE )

        #   add back the center of the zonohedron
        res = .Call( C_plusEqual, section, x$center, 1L )
        if( is.null(res) )  return(NULL)


        #   reorder the data
        section         = section[perm, ]
        idxhypersub     = idxhypersub[perm]
        antipodalmask   = antipodalmask[perm]


        df  = data.frame( row.names=rownames(section) )

        df$point    = section
        df$hyperidx = idxhypersub
        df$sign     = 1L - 2L * antipodalmask  #  so   TRUE -> -1L  and   FALSE -> +1L.  ifelse( antipodalmask, -1L, 1L )

        out[[k]]    = df    # list( beta=beta[k],  section=section )
        }

    if( FALSE )
        {
        timermain       = updatetimer( timermain )
        timesections    = timermain$elapsed
        timeall         = timermain$total

        cat( "edges:    ", timeedges * 1000, " msec\n" )
        cat( "facets:   ", timefacets * 1000, " msec\n" )
        cat( "centers:  ", timecenters * 1000, " msec\n" )
        cat( "sections: ", timesections * 1000, " msec",  " sections=", length(beta),  "   ", timesections*1000/length(beta), "msec per section\n" )
        cat( "total:    ", timeall * 1000, " msec\n" )
        }

    if( plot )
        {
        if( ! requireNamespace( 'rgl', quietly=TRUE ) )
            {
            log_level( WARN, "Package 'rgl' is required for plotting.  Please install it." )
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
                }
            }
        }

    return(out)
    }




#   given a point on the boundary, return a point in the unit n-cube that maps to it, under W
#       x               the zonohedron
#       point           Mx3 matrix with points on the boundary of x in the rows
#                       Such points typically come in 1 of 2 ways:
#                           1) as computed by raytrace()
#                           2) as computed by section()
#       tol             tolerance for being on the boundary
#
#   returns a data.frame with M rows and these columns:
#       point       the original given matrix of boundary points
#       distance    signed distance to the boundary of the zonohedron
#       facetidx    index of the facet pair
#       sign        sign of the facet in the pair
#       pcube       an MxN matrix, where M=nrow(boundarydata) and N=number of points in x$matroid
#                   each row of the matrix is a point in the n-cube
#                   that maps to the boundary point on the zonohedron
#       transitions the number of transitions - a non-negative even integer

invertboundary.zonohedron <- function( x, point, tol=5.e-14 )
    {
    point   = prepareNxM( point, 3 )

    if( is.null(point) )    return(NULL)

    ok  = is.numeric(point)  &&  is.matrix(point) &&  0<nrow(point)  &&  ncol(point)==3
    if( ! ok )
        {
        log_level( ERROR, "argument point is invalid." )
        return(NULL)
        }

    #   compute boundarydata from x and boundary
    direction   = .Call( C_sumMatVec, point, -x$center, 1L )

    #   in the next call, invert=FALSE
    boundarydata    = raytrace( x, x$center, direction )
    if( is.null(boundarydata) ) return(NULL)

    #   replace computed boundary points with the originals
    boundarydata$point   = point

    #print( boundarydata )

    #   check tmax
    #delta   = max( abs(boundarydata$tmax-1) )
    #if( tol < delta )
    #    {
    #    log.string( WARN, "boundary delta = %g > %g.", delta, tol )
    #    }

    dat = invertboundarydata( x, boundarydata, tol=tol )
    if( is.null(dat) )  return(NULL)

    out = data.frame( row.names=rownames(dat) )
    out$point       = point
    out$distance    = dat$distance
    out$facetidx    = boundarydata$facetidx
    out$sign        = boundarydata$sign
    out$pcube       = dat$pcube
    out$transitions = dat$transitions

    return( out )
    }


#   x               the zonohedron
#   boundarydata    as computed by raytrace()
#                   the only columns required are:  point, facetidx, sign
#
#   this function does the real work for invertboundary()
#   it is also called from raytrace(), when the arg invert=TRUE
#
#   returns a data.frame with the same number of rows, and new columns:
#       distance    signed distance to the boundary of the zonohedron
#       pcube       an MxN matrix, where M=nrow(boundarydata) and N=number of points in x$matroid
#                   each row of the matrix is a point in the n-cube
#                   that maps to the boundary point on the zonohedron
#       transitions the number of transitions - a non-negative even integer

invertboundarydata <- function( x, boundarydata, tol=5.e-14 )
    {
    ok  = is.data.frame(boundarydata)  &&  ! is.null(boundarydata$point)  &&  ! is.null(boundarydata$facetidx) &&  ! is.null(boundarydata$sign)
    if( ! ok )
        {
        log_level( ERROR, "data.frame argument boundarydata is invalid." )
        return(NULL)
        }

    point   = boundarydata$point

    matroidsimple   = getsimplified(x$matroid)
    matrixsimple    = getmatrix( matroidsimple )

    m       = nrow( boundarydata )
    nsimp   = ncol( matrixsimple )

    idxfromground   = idxfromgroundfun( matroidsimple$ground )

    distance    = rep( NA_real_, m )
    pcube       = matrix( NA_real_, m, nsimp )

    tolW    = 5.e-9
    tolE    = 5.e-5

    for( k in 1:m )
        {
        #   get facet center and normal, in the centered zonohedron
        facetidx    = boundarydata$facetidx[k]
        normal      = boundarydata$sign[k] * x$facet$normal[facetidx, ]
        center      = boundarydata$sign[k] * x$facet$center[facetidx, ]

        #   translate given point to centered zonohedron
        pcentered   = point[k, ] - x$center

        if( is.na(facetidx)  ||  is.na(pcentered[1]) )  next

        #   computed signed distance to the boundary
        distance[k] = sum( pcentered * normal )  -  x$facet$beta[ facetidx ]

        if( tol < abs(distance[k]) )    next    # too far from boundary


        colidx  = idxfromground[ matroidsimple$hyperplane[[ facetidx ]] ]
        #   cat( "colidx=", colidx, '\n' )


        numgen  = length( colidx )

        #   now compute the vector alpha, which has length numgen
        #   alpha has the coefficients of the edges of the facet

        if( numgen == 2 )
            {
            #   facet is a parallelogram, which is the usual case and easier case
            edge2   = matrixsimple[ , colidx ]     #; cat( "colidx=", colidx, '\n' )
            M       = cbind( edge2, normal )        #; print( M )    # M is 3x3
            y       = base::solve( M, pcentered - center )    #; print(y)
            y       = y[1:2]

            #   test for inside parallelogram, with a tolerance
            ok  = all( abs(y) <= 0.5 + tolW )
            if( ! ok )
                {
                if( all( abs(y) <= 0.5 + tolE ) )   # tolE is much bigger than tolW
                    {
                    lev = WARN      # keep going and translate and clamp
                    tol2 = tolW
                    }
                else
                    {
                    lev = ERROR     # this will force stoppage
                    tol2 = tolE
                    }

                log_level( lev, "Internal problem.  y=[%.15g,%.15g] is outside the square [-1/2,1/2]^2.  tol2=%g",
                                                    y[1], y[2], tol2 )
                # next    # something went wrong
                }

            # translate from [-0.5,0.5] to [0,1]
            alpha   = y + 0.5

            #   clamp alpha to [0,1]
            alpha   = pmin( pmax( alpha, 0), 1 )
            }
        else
            {
            #   facet is a more complex zonogon, with 6 or more sides
            zono2D  = x$zonogon[[ facetidx ]]

            #   check sizes
            if( numgen != length( getground( getsimplified(zono2D$matroid) ) ) )
                {
                log_level( ERROR, "Internal error.  numgen=%d  != %d.",
                                        numgen, length( getground( getsimplified(zono2D$matroid) ) ) )
                next    # something went wrong
                }

            A   = x$frame3x2[[ facetidx ]]
            p2D = as.double( (pcentered - center) %*% A )  +  zono2D$center

            # cat( "p2D=", p2D, '\n' )

            res = invert( zono2D, p2D, tol=tol )
            if( is.null(res) )   next   # something went wrong

            #   res$pcube has already been clamped to [0,1]

            alpha   = as.double( res$pcube )

            #   cat( "alpha=", alpha, '\n' )
            }


        #   assign the source point in the cube,
        #   except that the generators corresponding to the edges of this facet are wrong,
        #   and fixed by overriding with alpha in the next line
        pcube[ k, ]    = (sign(normal  %*%  matrixsimple) + 1) / 2

        #   override the coefficients of the edges of the facet
        pcube[ k, colidx ]   = alpha
        }

    #   print( pcube )

    pcube   = invertcubepoints( x, pcube, tol=tol )
    if( is.null(pcube) )    return(NULL)

    #   colnames(pcube) = as.character( x$matroid$ground )

    transitions = rep( NA_integer_, m )
    for( k in 1:m )
        {
        if( is.finite( pcube[k,1] ) )   transitions[k]  = transitioncount( pcube[k, ] )
        }

    rnames  = rownames(point)
    if( is.null(rnames) || anyDuplicated(rnames)!=0 )   rnames = 1:m

    out = data.frame( row.names=rnames )

    out$point       = point
    out$distance    = distance
    #out$facetidx    = boundarydata$facetidx
    #out$sign        = boundarydata$sign
    out$pcube       = pcube
    out$transitions = transitions

    if( FALSE )
        {
        #   test the inversion
        matorig = getmatrix( x$matroid )

        delta   = abs( pcube %*% t(matorig)  -  point)
        #cat( "range(delta)=", range(delta), '\n' )
        #print( t(matorig) )
        #print( delta )

        delta   = rowSums( delta )

        if( any( tol < delta, na.rm=TRUE ) )
            {
            log_level( WARN, "Inversion test failed.  max(delta)=%g > %g=tol",
                                max(delta,na.rm=TRUE), tol )
            }

        out$delta   = delta
        }

    return( out )
    }


#   x       a zonohedron object
#   gndpair Mx2 integer matrix, where rows define 2 points in the ground set of the matroid   #  *simplified* removed
#           these define 2 generators of the matroid, and a pgram in the boundary of the zonohedron.
#           If gndpair[1,j] < gndpair[2,j], the 'positive' facet is returned.
#           If gndpair[1,j] > gndpair[2,j], the 'negative' facet is returned.
#           If gndpair[1,j]==gndpair[2,j], the facet is undefined, and NA is returned
#   cube    if TRUE, then the corresponding point in the source cube is returned
#
#   returns a data.frame with M rows and these columns:
#       gndpair         the given gndpair
#       hyperplaneidx   the index of the hyperplane that contains gndpair
#       center          Mx3 matrix with center of the facet, relative to the center of the zonohedron
#       transitions     the number of transitions - a non-negative even integer
#   And if cube is TRUE, then this column
#       pcube       an MxN matrix, where M=nrow(boundarydata) and N=number of points in x$matroid$simplified
#                   each row of the matrix is a point in the n-cube
#                   that maps to the boundary point on the zonohedron
#
#   the row.names are set to the index of the corresponding gndpair

boundarypgramdata <- function( x, gndpair, cube=FALSE )
    {
    if( ! inherits(x,"zonohedron") )
        {
        log_level( ERROR, "Argument x is invalid.  It is not a zonohedron." )
        return(NULL)
        }

    gndpair = prepareNxM( gndpair, 2 )
    if( is.null(gndpair) )  return(NULL)

    dimsave     = dim(gndpair)

    if( ! is.integer(gndpair) )
        {
        #   it's probably floating point, because the user typed gndpair on the command line
        #   change to true integer for the output
        gndpair = as.integer(gndpair)
        dim(gndpair) = dimsave          # now it can be assigned to the returned data.frame
        }


    ground  = getmatroid(x)$ground      # ground set of original matroid

    #   change any entries in gndpair, that are not in ground, to NA_integer_
    gndpairsave = gndpair               # because we might modify gndpair with NAs
    mask    = ! (gndpair %in% ground)
    if( any(mask) )
        {
        gndpair[mask]   = NA_integer_
        dim(gndpair)    = dimsave
        }

    #   convert gndpair to raw indexes in the original matroid
    idxfromground   = idxfromgroundfun( ground )      # idxfromgroundfun( matsimple$ground )

    idxpair_org         = idxfromground[ gndpair ]
    dim(idxpair_org)    = dim(gndpair)


    #   collapse to raw indexes in the simplified matroid
    collapsetosimple    = getmatroid(x)$collapsetosimple

    if( is.null(collapsetosimple) ) collapsetosimple = 1:length(ground)   # matroid is already simple

    idxpair = collapsetosimple[ idxpair_org ]
    dim( idxpair )  = dim( gndpair )        # ; print( idxpair )

    #   idxpair now contains raw indexes in the simplified matroid
    #   because of the collapse, there might be duplicates in the rows, even if there were not there before

    #   swap pairs if necessary
    signvec = sign( idxpair[ ,2] - idxpair[ ,1] )       # signvec might contain NAs
    mask    = signvec < 0
    mask[ is.na(mask) ] = FALSE
    if( any(mask) )
        {
        idxpair[mask, ]     = idxpair[mask, 2:1]
        idxpair_org[mask, ] = idxpair_org[mask, 2:1]
        }


    matsimple   = getsimplified( getmatroid(x) )

    idxfromgroundsimple = idxfromgroundfun( matsimple$ground )      # idxfromgroundfun( matsimple$ground )

    m       = nrow( idxpair )

    n       = length( ground )
    nsimp   = length( matsimple$ground )

    #   in the next line, invalid idxpair entries lead to NA values in pairidx, which is a 1-based integer vector of length m
    #   in particular, if idxpair[ ,1] == idxpair[ ,2] then the value of pairidx is NA_integer_
    pairidx = .Call( C_pairindex, idxpair, nsimp )  #; print( pairidx )

    matrix          = getmatrix( getmatroid(x) )

    #   matrixsimple    = getmatrix( matsimple )

    loopindexes     = idxfromground[ getmatroid(x)$loop ]

    #   NA values in pairidx lead to NA values in hyperplaneidx
    hyperplaneidx   = matsimple$hyperplaneidx[ pairidx ]

    pcube           = matrix( NA_real_, m, n )

    transitions     = rep( NA_integer_, m )

    for( k in 1:m )
        {
        # cat( "------------  k=", k, '---------\n' )

        hyperidx    = hyperplaneidx[k]      # = matsimple$hyperplaneidx[ pairidx[k] ]

        if( is.na(hyperidx) )   next

        normal  = x$facet$normal[hyperidx, ]

        pc  = normal  %*%  matrix       #;       print( pc[colidx] )


        #   get the raw indexes in the simplified matroid for this facet
        colidx  = idxfromgroundsimple[ matsimple$hyperplane[[ hyperidx ]] ]  #; print( colidx )

        #print( "colidx" )
        #print( colidx )

        #   lift the indexes to raw index in the original matrix
        colidx  = liftrawindexes( getmatroid(x), colidx )

        #   add any loops to colidx, the pc[loop] should be exactly 0
        colidx  = c( colidx, loopindexes )


        #   pc[colidx] should be 0 or nearly 0.
        #   override and force to exactly 0
        pc[colidx] = 0

        if( FALSE )
            {
            #   check that the zeros of pc[] are exactly colidx
            idx = which(pc==0)
            ok  = length(idx)==length(colidx)  &&  all( idx == colidx )
            if( ! ok )
                {
                log_level( WARN, "internal error.  colidx  !=  which(pc==0)" )

                print( "colidx" )
                print( colidx )

                print( "which(pc==0)" )
                print( which(pc==0) )
                }
            }

        #   scale sign from [-1,+1] to [0,1]
        pc  = (sign(pc) + 1) / 2        #;   print( pc )

        if( 2 < length(colidx) )
            {
            #   the generators are part of a non-pgram zonogon face
            #   make changes to colidx entries
            changeable  = logical(n)
            changeable[ colidx ]    = TRUE
            changeable[ idxpair_org[k, ] ]  = FALSE     # but do not change the 2 given indexes

            #   set "inside" to 1, and "outside" to 0
            inside  = logical(n)
            inside[ idxpair_org[k,1]:idxpair_org[k,2] ] = TRUE

            pc[  inside & changeable ]  = 1
            pc[ !inside & changeable ]  = 0

            #   now only only 2 entries of pc[] are equal to 1/2, which is what we want for a parallelogram
            }

        transitions[k]  = transitioncount( pc )

        if( isTRUE(0 < signvec[k]) )
            pcube[k, ] = pc
        else
            pcube[k, ] = 1 - pc     # complement
        }


    rnames  = rownames( gndpairsave )
    if( is.null(rnames) )   rnames  = 1:m

    out = data.frame( row.names=rnames )

    out$gndpair         = gndpairsave
    out$hyperplaneidx   = hyperplaneidx
    out$center          = tcrossprod( pcube, matrix )   # pcube  %*%  t(matrix)      #         x$facet$center[ hyperplaneidx, , drop=FALSE]
    out$transitions     = transitions

    if( cube )  out$pcube   = pcube     # add cube points


    #   fix the signs of the center
    #   out$center  = signvec * out$center  #   all columns are multiplied by signvec

    if( FALSE  &&  cube )
        {
        #   'lift' cube points from the simplified matroid dimension to the original matroid dimension
        out$pcube   = invertcubepoints( x, pcube )

        for( k in 1:m )
            {
            if( signvec[k] < 0 )
                #   invert the k'th point
                out$pcube[k, ] = 1 - out$pcube[k, ]
            }
        }

    return( out )
    }






if( FALSE )
{
#   x           a zonohedron object
#   direction   Mx3 matrix, with the M directions in the rows, direction (0,0,0) is invalid
#   tol         tolerance for argmax, being in the same affine subspace
#
#   returns a data.frame with M rows and these columns:
#       direction   the given matrix of directions
#       value       the value of the support function of x, in the given direction
#       argmax      a point on the boundary of x where the max is taken
#       dimension   of the set argmax, 0 means a vertex and 1 means a facet
#
#   value:  see support_zonotope()

support.zonohedron <- function( x, direction, tol=5.e-15 )
    {
    return( support_zonotope(x,direction,tol) )
    }

#   x   a zonohedron whose matroid is simple
#
#   replace each generator g by the pair -g/2 , +g/2

symmetrize.zonohedron <- function( x )
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

    out = zonohedron( matgen, e0=0, e1=0, ground=gndgen )

    return(out)
    }


minkowskisum.zonohedron  <- function( zono1, zono2, e0=0, e1=1.e-6, e2=1.e-10, ground=NULL, ... )
    {
    if( ! inherits(zono2,"zonohedron") )
        {
        log_level( ERROR, "2nd argument zono2 is invalid.  It's not a zonohedron." )
        return(NULL)
        }

    #   get the 2 matrices and cbind them
    mat1    = zono1$matroid$matrix
    mat2    = zono2$matroid$matrix

    out = zonohedron( cbind(mat1,mat2), e0=e0, e1=e1, e2=e2, ground=ground )

    return( out )
    }

'%+%.zonohedron'  <- function(zono1,zono2)
    {
    return( minkowskisum( zono1, zono2 ) )
    }
}

if( FALSE )
{
is_pointed.zonohedron <- function( x )
    {
    return( 3 <= length(x$facet0) )
    }

is_salient.zonohedron <- function( x )
    {
    return( 0 < length(x$facet0) )
    }
}



#   x       a zonohedron
#   gen     point in the ground set refering to a generator
#   full    if TRUE, return both halves of the belt, if FALSE then only return half the belt
#   returns a dataframe with 2 columns
#       *) point0   m x 2 matrix of edge start points
#       *) point1   m x 2 matrix of edge stop points

getbeltedges <- function( x, gen, full=TRUE )
    {
    matsimp = getsimplified( x$matroid )

    gen = as.integer(gen)  #;  print( str(gen) )

    genidx  = match( gen, matsimp$ground )

    #   midpointmat is for the centered zonohedron
    midpointmat = .Call( C_beltmidpoints, matsimp$hyperplane, x$beltlist[[genidx]], gen, x$facet$center, x$facet$normal,
                                            matsimp$ground, matsimp$matrix, matsimp$crossprods )
    if( is.null(midpointmat) )    return(NULL)

    if( full )  midpointmat = rbind( midpointmat, -midpointmat )

    .Call( C_plusEqual, midpointmat, x$center, 1L )

    edge    = matsimp$matrix[ ,genidx]

    #point0  = duplicate(midpointmat)
    #.Call( C_plusEqual, point0, -0.5*edge, 1L )
    point0  = .Call( C_sumMatVec, midpointmat, -0.5*edge, 1L )

    #point1  = duplicate(midpointmat)
    #.Call( C_plusEqual, point1,  0.5*edge, 1L )
    point1  = .Call( C_sumMatVec, midpointmat, 0.5*edge, 1L )


    out = data.frame( row.names=1:nrow(midpointmat) )
    out$midpointmat = midpointmat
    out$point0      = point0
    out$point1      = point1

    return(out)
    }

if( FALSE )
{
#   methods taken from:
#       Optimal Whitening and Decorrelation
#       Agnan Kessy1, Alex Lewin, and Korbinian Strimmer (2016)
#
#   returns a new zonohedron, as spherical as possible

spherize.zonohedron <- function( x, method="ZCA", ... )
    {
    return( spherize_zonotope( x, method=method ) )
    }
}



#   x   a zonohedron
#   W   a 3x3 invertible matrix.  Put matrix on the left, and vector on the right.

lintransform.zonohedron <- function( x, W )
    {
    if( length(W) == 1 )
        W = W * diag(3)

    ok  = all( dim(W) == c(3,3) )
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

    if( is.null(out$matroid) )  return(NULL)

    out$center  =   as.double( x$center %*% t(W) )

    out$facet$center    = x$facet$center %*% t(W)

    #   compute all the outward pointing edge normals
    normal  = x$facet$normal %*% Winv

    #   now unitize
    normal  = normalizeMatrix( normal, 1L )
    out$facet$normal  = normal

    #   calculate the n plane constants beta
    #   these plane constants are for the centered zonogon
    out$facet$beta    = .rowSums( normal * out$facet$center, nrow(normal), ncol(normal) )

    if( ! is.null(out$frame3x2) )
        {
        for( k in 1:length(out$frame3x2) )
            out$frame3x2[[k]] = W %*% out$frame3x2[[k]]
        }

    if( ! is.null(out$zonoseg) )
        {
        #   recompute out$zonoseg list completely
        out$zonoseg     = makezonoseglist( out$matroid )
        }

    attr( out, "lintransform" ) = W

    return( out )
    }


initplot3D <- function( zono, bgcol )
    {
    center  = zono$center
    white   = 2 * center

    #   start 3D drawing
    rgl::bg3d( color=bgcol )

    rgl::points3d( 0, 0, 0, col='black', size=10, point_antialias=TRUE )
    rgl::points3d( white[1], white[2], white[3], col='white', size=10, point_antialias=TRUE )
    rgl::points3d( center[1], center[2], center[3], col='gray50', size=10, point_antialias=TRUE )

    #   exact diagonal of zono
    rgl::lines3d( c(0,white[1]), c(0,white[2]), c(0,white[3]), col=c('black','white'), lwd=3, lit=FALSE )

    return(TRUE)
    }


#   x       a zonohedron object
#   type    'e' for edges, 'p' for points drawn at the center of each p-face, 'f' for filled faces
#   both    draw both symmetric halves


plot.zonohedron <- function( x, type='e', pcol=NULL, ecol=NULL, ewd=3, etcol=NA,
                                fcol=NULL, falpha=1, normals=FALSE, bgcol="gray40", both=TRUE, ... )
    {
    if( ! requireNamespace( 'rgl', quietly=TRUE ) )
        {
        log_level( ERROR, "Package 'rgl' cannot be loaded. It is required for plotting the zonohedron." )
        return(FALSE)
        }

    center  = x$center
    white   = 2 * center


    #   start 3D drawing
    rgl::bg3d( color=bgcol )


    cube    = rgl::scale3d( rgl::cube3d(col="white"), center[1], center[2], center[3] )
    cube    = rgl::translate3d( cube, center[1], center[2], center[3]  )


    rgl::points3d( 0, 0, 0, col='black', size=10, point_antialias=TRUE )
    rgl::points3d( white[1], white[2], white[3], col='white', size=10, point_antialias=TRUE )
    rgl::points3d( center[1], center[2], center[3], col='gray50', size=10, point_antialias=TRUE )

    #   exact diagonal of box
    rgl::lines3d( c(0,white[1]), c(0,white[2]), c(0,white[3]), col=c('black','white'), lwd=3, lit=FALSE )

    matsimp = getsimplified( x$matroid )

    matgen  = getmatrix(matsimp)
    numgen  = ncol(matgen)

    gndgen  = getground(matsimp)

    edgecoeff   = matrix( c( -0.5,-0.5, -0.5,0.5, 0.5,0.5, 0.5,-0.5), 2, 4 )    # 2x4

    if( grepl( 'e', type ) )
        {
        #   wireframe
        # rgl::wire3d( cube, lit=FALSE )

        n   = length(gndgen)

        if( is.null(ecol) )
            colvec  = rainbow( n )
        else
            {
            colvec  = ecol
            m   = length(colvec)
            if( m < n ) colvec  = c( colvec, rep(colvec[m],n-m) )   #   extend with the last color
            }

        for( i in 1:n )
            {
            edgedf  = getbeltedges( x, gndgen[i], full=both )

            xyz = rbind( edgedf$point0, edgedf$point1 )

            m   = nrow(edgedf)

            perm    = 1:(2*m)
            dim(perm)    = c(m,2L)
            perm = t(perm)
            dim(perm)    = NULL
            # print( perm )

            xyz = xyz[ perm, ]

            rgl::segments3d( xyz, col=colvec[i], lwd=ewd )    #, front=polymode, back=polymode, col=col, lit=FALSE )
            }

        facetsNT    = length(x$zonogon)     # the number of Non-Trivial facets in the zonohedron

        if( ! is.na(etcol)  &&  0<facetsNT )
            {
            #   draw the pgram edges for the tiling of each non-pgram facet

            for( i in seq_len(facetsNT) )
                {
                zono    = x$zonogon[[i]]

                matrixzono  = getmatrix( getsimplified(zono$matroid) )

                #   rotate into 3D
                matrixzono = x$frame3x2[[i]]  %*%  matrixzono

                #   find center of all pgram tiles in 3D
                center3D    = zono$tilingdata$center %*% t( x$frame3x2[[i]] )

                #   translate each pgram center to the zonogon facet in 3D
                center3D    = .Call( C_sumMatVec, center3D, x$facet$center[i, ], 1L )

                pgrams  = nrow(center3D)
                step    = 4
                quadmat = matrix( 0, nrow=step*pgrams, ncol=3 )

                for( j in 1:pgrams )
                    {
                    col2    = zono$tilingdata$idxpair[j, ]

                    edge    = matrixzono[ , col2 ]    # 3x2

                    k       = step*(j-1)

                    #quadmat[k+1, ] = center3D[j, ]  - 0.5 * edge[ , 1] - 0.5*edge[ , 2]
                    #quadmat[k+2, ] = center3D[j, ]  - 0.5 * edge[ , 1] + 0.5*edge[ , 2]
                    #quadmat[k+3, ] = center3D[j, ]  + 0.5 * edge[ , 1] + 0.5*edge[ , 2]
                    #quadmat[k+4, ] = center3D[j, ]  + 0.5 * edge[ , 1] - 0.5*edge[ , 2]

                    quadmat[ (k+1):(k+4),  ] =  .Call( C_sumMatVec, t(edge %*% edgecoeff), center3D[j, ], 1L )
                    }

                xyz = .Call( C_sumMatVec, quadmat, x$center, 1L )

                if( both )
                    xyz = rbind( xyz, .Call( C_sumMatVec, -quadmat, x$center, 1L ) )

                rgl::quads3d( xyz, col=etcol, lwd=1, front='lines', back='lines', lit=FALSE  )   # quad edges
                }
            }
        }


    if( grepl( 'p', type ) )
        {
        colvec  = pcol

        if( is.null(colvec) )
            colvec  = c( 'black', 'red' )
        else if( length(colvec) == 1 )
            colvec  = rep( colvec[1], 2 )

        #   draw first half in 'black'
        #xyz = duplicate( x$facet$center )
        #.Call( C_plusEqual, xyz, x$center, 1L )
        xyz = .Call( C_sumMatVec, x$facet$center, x$center, 1L )

        rgl::points3d( xyz[ ,1],  xyz[ ,2], xyz[ ,3], col=colvec[1], size=6, point_antialias=TRUE )

        if( both )
            {
            #   draw 2nd half in 'red'
            #xyz = duplicate( -(x$facet$center) )
            #.Call( C_plusEqual, xyz, x$center, 1L )
            xyz = .Call( C_sumMatVec, -(x$facet$center), x$center, 1L )

            rgl::points3d( xyz[ ,1],  xyz[ ,2], xyz[ ,3], col=colvec[2], size=6, point_antialias=TRUE )
            }
        }


    if( grepl( 'f', type ) )
        {
        #   draw filled quads
        #   make lookup table from ground to column index
        idxfromground   = integer( gndgen[numgen] )
        idxfromground[ gndgen ] = 1:numgen

        hyper   = matsimp$hyperplane
        lenvec  = lengths(hyper)

        idx2    = which( lenvec==2 )
        pgrams  = length(idx2)

        step    = 4
        quadmat = matrix( 0, nrow=step*pgrams, ncol=3 )

        if( is.null(fcol) )
            fcol    = c( 'blue', 'red', 'yellow', 'green', 'orange', 'purple' )

        #   plot the parallelograms

        for( j in 1:pgrams )
            {
            i   = idx2[j]   # index into hyper *and* facet

            center  = x$facet$center[i, ]

            edge    = matgen[  , idxfromground[ hyper[[i]] ] ]  # 3x2 matrix

            k       = step*(j-1)

            #quadmat[k+1, ] = center - 0.5 * edge[ , 1] - 0.5*edge[ , 2]
            #quadmat[k+2, ] = center - 0.5 * edge[ , 1] + 0.5*edge[ , 2]
            #quadmat[k+3, ] = center + 0.5 * edge[ , 1] + 0.5*edge[ , 2]
            #quadmat[k+4, ] = center + 0.5 * edge[ , 1] - 0.5*edge[ , 2]

            quadmat[ (k+1):(k+4),  ] =  .Call( C_sumMatVec, t(edge %*% edgecoeff), center, 1L )
            }


        xyz = .Call( C_sumMatVec, quadmat, x$center, 1L )

        if( both )
            #   add opposite half
            xyz = rbind( xyz, .Call( C_sumMatVec, -quadmat, x$center, 1L ) )

        rgl::quads3d( xyz, col=fcol[1], alpha=falpha )    #, front=polymode, back=polymode, col=col, lit=FALSE )



        #   plot the non-parallelograms

        for( k in seq_len( length(x$zonogon) ) )
            {
            zonok   = x$zonogon[[k]]

            #   center the vertices of the zonogon
            vertex  = .Call( C_sumMatVec, zonok$vertex, -zonok$center, 1L )

            #   map from 2D to 3D, center remains at 0
            vertex  = vertex %*% t(x$frame3x2[[k]])

            #   add the facet center and zonohedron center
            xyz = .Call( C_sumMatVec, vertex, x$center + x$facet$center[k, ], 1L )

            jmax    = which.max( abs(x$facet$normal[k, ]) )
            coord   = 1:3
            coord   = coord[-jmax]

            numgen  = ncol( zonok$matroid$matrix )
            col     = fcol[ min( numgen-1, length(fcol) ) ]

            quadmat = makequads( xyz )

            rgl::quads3d( quadmat, col=col, alpha=falpha )

            if( both )
                {
                #   draw opposite half
                xyz = .Call( C_sumMatVec, vertex, x$center - x$facet$center[k, ], 1L )
                quadmat = makequads( xyz )
                rgl::quads3d( quadmat, col=col, alpha=falpha )

                #rgl::polygon3d( xyz[ ,1], xyz[ ,2], xyz[ ,3], fill=TRUE, coord=coord, col=col, random=FALSE )
                }
            }
        }


    if( normals )
        {
        xyz = .Call( C_sumMatVec, x$facet$center, x$center, 1L )

        for( i in 1:nrow(xyz) )
            rgl::arrow3d( xyz[i, ], xyz[i, ] + x$facet$normal[i, ], type="lines", col="black" )

        if( both )
            {
            xyz = .Call( C_sumMatVec, -x$facet$center, x$center, 1L )

            for( i in 1:nrow(xyz) )
                rgl::arrow3d( xyz[i, ], xyz[i, ] - x$facet$normal[i, ], type="lines", col="black" )
            }
        }

    return( invisible(TRUE) )
    }


as.mesh3d.zonohedron  <-  function( x, fcolor=NULL, falpha=1, codes=FALSE, ... )
    {
    #timermain = createtimer()

    if( is.null(fcolor) )   fcolor=c('blue','red','yellow','green', 'orange', 'purple')

    matsimp = getsimplified( x$matroid )

    matgen  = getmatrix(matsimp)
    numgen  = ncol(matgen)

    gndgen  = getground(matsimp)

    #   make lookup table from ground to column index
    idxfromground   = integer( gndgen[numgen] )
    idxfromground[ gndgen ] = 1:numgen

    hyper   = matsimp$hyperplane
    lenvec  = lengths(hyper)

    idx2    = which( lenvec==2 )
    pgrams  = length(idx2)

    facets      = 2L * nrow( x$facet )  # doubled, because we do the entire surface, not just half

    edges       = 2L * sum( lenvec )    # doubled, because we do the entire surface, not just half

    vertices    = edges - facets + 2L   # total vertices on boundary of x

    #   make C++ map to hold the cube vertex points, aka the "binary codes"
    handle  = .Call( C_makeRawMap, numgen, vertices )
    if( is.null(handle) )   return(NULL)

    on.exit( .Call( C_deleteRawMap, handle ) )      # this *does* work

    #   allocate matrix to hold all the vertices on boundary of x
    #  vertex  = matrix( NA_real_, nrow=3, ncol=vertices )

    #   allocate matrix to hold the integer quad indexes, into vertex[]
    #   quadpram is only for the pgram facets
    quadpgram   = matrix( NA_integer_, nrow=4, ncol=2L*pgrams )

    #cat( "-------------------\n" )
    #cat( "starting &vertex:", obj_addr(vertex), "   starting quadpgram:", obj_addr(quadpgram), '\n' )


    colorvec    = rep( fcolor[1], ncol(quadpgram) )

    #   make 2x4 matrices for the 4 pgram vertices
    edgecoeff   = matrix( c( -0.5,-0.5, -0.5,0.5, 0.5,0.5, 0.5,-0.5), 2, 4 )    # 2x4 matrix of coefficients, relative to the facet center

    #   make matrix for the 4 pgram vertices
    raw4        = as.raw( 0 < edgecoeff )           # + goes to raw 01, and - goes to raw 00
    dim(raw4)   = dim(edgecoeff)                    # 2x4

    # idx = integer(4)    # indexes for a single quad

    #timermain   = updatetimer( timermain )
    #timeprep    = timermain$elapsed


    #   process the parallelograms

    #timeindex   = 0

    vidxmax = 0     # max vertex index so far

    log_level( TRACE, "Processing %d parallelogram facets, and creating %d quads...", 2L*pgrams, ncol(quadpgram) )

    for( j in 1:pgrams )
        {
        i   = idx2[j]   # index into both hyper *and* facet

        # center      = x$facet$center[i, ]       #  center of the facet

        pgramgen    = idxfromground[ hyper[[i]] ]         #   the 2 column indexes for the pgram generators

        # edge        = matgen[  , pgramgen ]     # 3x2 matrix = the 2 edge generators of the pgram

        normal      = x$facet$normal[i, ]       # outward pointing unit normal of facet i

        sign        = x$facet$sign[i]

        #   make raw matrix of 0s and 1s, from which we will make 4 cube vertices
        #   the columns of pcuberaw are replicated
        pcuberaw    = array( as.raw( 0 < (normal %*% matgen) ), dim=c(numgen,4) )

        #   the 2 values pcuberaw[ pgramgen, ] are indeterminate,
        #   because they are very close to 0, but not exactly because there is floating point truncation
        #   but we now fill in these indeterminates, by assigning these indexes in all 4 columns, i.e. the 4 vertices of the pgram

        pcuberaw[ pgramgen, ]   = raw4

        #thetime = microbenchmark::get_nanotime()

        idx = .Call( C_getIndexRaw, handle, pcuberaw, FALSE )   # idx has length 4

        if( is.null(idx) )  return(NULL)

        #timeindex   = timeindex + microbenchmark::get_nanotime() - thetime

        if( any( vertices < idx ) )
            {
            log_level( FATAL, "Internal error.  idx = %d,%d,%d,%d > %d = number of vertices.",
                            idx[1], idx[2], idx[3], idx[4], vertices )
            return( NULL )
            }

        if( 0 < sign )  idx =  rev(idx)     # reverse order

        quadpgram[ , j] = idx



        #   repeat, but do the antipodal facets

        #thetime =  microbenchmark::get_nanotime()

        idx_anti    = .Call( C_getIndexRaw, handle, pcuberaw, TRUE )    # TRUE means take complement

        if( is.null(idx_anti) )  return(NULL)

        #timeindex   = timeindex + microbenchmark::get_nanotime() - thetime

        if( any( vertices < idx_anti ) )
            {
            log_level( FATAL, "Internal error.  idx_anti = %d,%d,%d,%d > %d = number of vertices.",
                            idx_anti[1], idx_anti[2], idx_anti[3], idx_anti[4], vertices )
            return( NULL )
            }

        if( sign < 0 )  idx_anti = rev(idx_anti)     # reverse order

        quadpgram[ , j+pgrams]  = idx_anti

        #  cat( "j=", j, "  sign=", sign,  "   idx=", idx, "\n" )

        vidxmax = max( vidxmax, idx, idx_anti )
        }

    log_level( TRACE, "... added %d of %d vertices.", vidxmax, vertices )

    #cat( "ending   &vertex:", obj_addr(vertex), "   ending   quadpgram:", obj_addr(quadpgram), '\n' )

    #cat( "getindex:  ", 1.e-9 * timeindex, " sec\n" )


    #timermain   = updatetimer( timermain )
    #timepgrams  = timermain$elapsed


    #   process the non-parallelograms, with 3 or more generators
    zonogons    = length(x$zonogon)

    quads   = 2L * ifelse( 0 < zonogons, sum( lenvec[1:zonogons] ) - zonogons, 0L )

    quad    = matrix( NA_integer_, nrow=4, ncol=quads )

    quadcount   = 0L

    # idx = integer(4)    # indexes for a single quad

    log_level( TRACE, "Processing %d non-trivial zonogon facets, and creating %d quads ...", 2L*zonogons, quads )

    for( j in seq_len(zonogons) )
        {
        #cat( "--------------  zonogon ", j, " of ", zonogons, '  -----------\n' )

        zonoj   = x$zonogon[[j]]

        # center  = x$facet$center[j, ]       #  center of the facet

        # sign        = x$facet$sign[j]


        #   let n be the number of generators of this facet.  n >= 3
        #   the number of vertices in this facet is 2*n
        n   = ncol( zonoj$matroid$matrix )

        if( n != length(hyper[[j]]) )
            {
            log_level( FATAL, "n = %d  !=  %d = length(hyper[[%d]]).", n, length(hyper[[j]]), j )
            return(NULL)
            }

        normal      = x$facet$normal[j, ]       # outward pointing unit normal of facet i

        zonojgen    = idxfromground[ hyper[[j]] ]   #   the n column indexes in matgen[,] for the generators of zonoj

        # edge        = matgen[  , zonojgen ]     # 3xN matrix = the N edge generators of the zonogon

        #   make numgen x 2n raw matrix of 0s and 1s, from which we will make the vertices of the facet
        pcuberaw    = array( as.raw( 0 < (normal %*% matgen) ), dim=c(numgen,2L*n) )    # columns are replicated  2*n times

        #   change entries in the generator rows - zonojgen - which will make all the columns of pcuberaw distinct
        pcuberaw[ zonojgen, ] = t( zonoj$pcube )      #  n x 2n


        #   get the indexes of the 2n vertices in the facet
        idx = .Call( C_getIndexRaw, handle, pcuberaw, FALSE )   # FALSE means do not do the antipodal facet

        #   idx should have length 2n
        if( length(idx) != 2L*n )
            {
            log_level( FATAL, "Internal error.  length(idx)=%d, but expected %d.", length(idx), 2L*n );
            return(NULL)
            }

        #   convert 0-1 raw to +/- 0.5 floating-point
        #edgecoeff       = as.double( pcube ) - 0.5
        #dim(edgecoeff)  = dim( pcube )   #  N x 2N,

        #   tile this zonogon with N-1 quads
        quadidx = makequadindexes( n )     # (N-1) x 4


        for( i in 1:nrow(quadidx) )
            {
            #   add column to quad[]
            quadcount   = quadcount + 1L

            if( ncol(quad) < quadcount )
                {
                log_level( FATAL, "Internal error.  quad index = %d > %d = allocated number of quads.",
                                quadcount, ncol(quad) )
                return( NULL )
                }

            quad[ , quadcount]  = idx[ quadidx[i, ] ]      # do not reverse here
            }


        #   repeat, but this time do the antipodal quads

        #   get the indexes of the 2n vertices in the antipodal facet
        idx_anti    = .Call( C_getIndexRaw, handle, pcuberaw, TRUE )   # TRUE means do the antipodal facet

        for( i in 1:nrow(quadidx) )
            {
            #   add column to quad[]
            quadcount   = quadcount + 1L

            if( ncol(quad) < quadcount )
                {
                log_level( FATAL, "Internal error.  quad index = %d > %d = allocated number of quads.",
                                quadcount, ncol(quad) )
                return( NULL )
                }

            quad[ , quadcount]  = rev( idx_anti[ quadidx[i, ] ] )      # reverse order for antipodals
            }

        vidxmax = max( vidxmax, idx, idx_anti )

        #   and finally, append the colors for these quads
        k   = min( n-1, length(fcolor) )

        colorvec    = c( colorvec, rep(fcolor[k],2*nrow(quadidx)) ) # 2 because we also make the antipodal facet
        }

    log_level( TRACE, "Added %d of %d vertices, and %d quads.", vidxmax, vertices, quadcount )

    if( vidxmax != vertices )
        {
        log_level( FATAL, "Internal error.  Created %d vertices, but expected to create %d.",
                                vidxmax, vertices )
        return(NULL)
        }

    if( quadcount != ncol(quad) )
        {
        log_level( FATAL, "Internal error.  Created %d quads, but expected to create %d.",
                                quadcount, ncol(quad) )
        return(NULL)
        }

    #timermain   = updatetimer( timermain )
    #timezonos   = timermain$elapsed



    #  .Call( C_deleteRawMap, handle )  handled OK in on.exit()

    #   vertex[] is now centered.  translate all vertices by x$center in place
    # .Call( C_plusEqual, vertex, x$center, 2L )

    #   the map, pointed to by handle, is now full of binary codes
    #   compute all of vertices from these codes, no multiplication is required, only addition
    vertex  = .Call( C_computeVertices, handle, matgen )


    out = list()

    out$vb          = rbind(vertex,1)                       #   for vb, add a row of 1s to make homogeneous coords
    out$ib          = cbind(quadpgram,quad)                 #   for ib, cbind the pgram and the non-pgram quads
    out$material    = list( color=colorvec, alpha=falpha )
    out$meshColor   = 'faces'

    if( codes ) out$codes   = .Call( C_getCodes, handle )

    class(out)  = c( "mesh3d", "shape3d" )


    if( FALSE )
        {
        timermain   = updatetimer( timermain )
        timefinish  = timermain$elapsed

        timetotal   = timermain$total

        cat( "preparation:  ", timeprep, " sec\n" )
        cat( "pgrams        ", timepgrams, " sec\n" )
        cat( "zonogons:     ", timezonos, " sec\n" )
        cat( "finish        ", timefinish, " sec\n" )
        cat( "total:        ", timetotal, " sec\n" )
        }



    return( out )
    }




plotpolygon <- function( x, normal=NULL, points=TRUE, labels=TRUE )
    {
    if( ! inherits(x,"zonohedron") )
        {
        log_level( ERROR, "Argument x is invalid.  It's not a zonohedron." )
        return(NULL)
        }

    if( ! is_pointed(x) )
        {
        log_level( ERROR, "Cannot plot polygon because the zonohedron is not pointed." )
        return( FALSE )
        }


    if( is.null(normal) )
        {
        #   make a matrix of candidates to try
        normalmat = matrix( c(1,1,1, 0,0,1, 0,1,0, 1,0,0), ncol=3, byrow=TRUE )
        normalmat = rbind( normalmat, supportingnormal0(x) )
        }
    else
        {
        ok  = is.numeric(normal)  &&  length(normal)==3  &&  any( normal!=0 )
        if( ! ok )
            {
            log_level( ERROR, "Argument normal is invalid." )
            return(FALSE)
            }

        normalmat = matrix( normal, nrow=1 )
        }

    genmat  = getmatrix( getsimplified(x$matroid) )  # 3 x N

    #   try all rows in normalmat
    found   = FALSE
    for( i in 1:nrow(normalmat) )
        {
        dots    = normalmat[i, ,drop=F] %*% genmat
        #   print( i ) ; print( dots )

        if( all( 0 < dots ) )
            {
            #   verified
            found   = TRUE
            normal  = normalmat[i, ]
            break
            }
        }

    if( ! found )
        {
        if( is.null(normal) )
            log_level( ERROR, "Cannot find a valid normal vector and halfspace." )
        else
            log_level( ERROR, "The given normal=%g,%g,%g is invalid.", normal[1], normal[2], normal[3] )
        return( FALSE )
        }

    #   transpose and rescale genmat
    genmat = t(genmat) / as.double(dots)       # N x 3  using recycling rule

    if( all(normal==c(1,1,1)) )
        frame3x2    = matrix( c(1,0,0, 0,1,0), nrow=3, ncol=2 )
    else
        frame3x2 = frame3x2fun( normal, TRUE )

    uv  = genmat %*% frame3x2   # N x 2

    xlim    = range( uv[ ,1] )
    ylim    = range( uv[ ,2] )
    plot( xlim, ylim, type='n', asp=1, lab=c(10,10,7), xlab='u', ylab='v', las=1 )
    grid( lty=1 )
    abline( h=0, v=0 )

    polygon( uv[ ,1], uv[ ,2], col=NA, border='black' )

    if( points )    points( uv[ ,1], uv[ ,2], pch=20, cex=0.8 )

    if( labels )
        {
        text( 1.05 * uv[ ,1], 1.05 * uv[ ,2], rownames(genmat) )
        }

    zononame = deparse(substitute(x))

    main = sprintf( "%s\n normal=[%g,%g,%g]", zononame, normal[1], normal[2], normal[3] )
    title( main=main, cex.main=0.8 )

    return( invisible(TRUE) )
    }




print.genlist <- function( x, full=TRUE, ... )
    {
    zlist   = lapply( x, zonohedron )

    # names(zlist)    = NULL
    #print( do.call( summary, zlist ) )

    print( summary_from_zlist( zlist, full=full ) )

    return( invisible(TRUE) )
    }

