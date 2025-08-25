
#   matroid is an S3 class and a list with these items:
#
#   ground          a positive integer vector in ascending order giving the ground set
#   hyperplane      a list of integer vectors - the hyperplanes that define the matroid.  All are subsets of the ground set.
#   rank            an integer - the rank
#   loop            an integer vector - the loops - each is a point in the ground set.  Can be empty.
#   multiple        a list of integer vectors - each is a group of (non-trivial) multiple points in ground set.  Can be empty.
#   matrix          the real matrix of generators - only if the matroid is constructed from a matrix
#   multiplesupp    a data.frame with #rows = length(multiple), and these columns:
#           contiguous  are the points of the group contiguous in the ground set minus loops, in the cyclic sense
#       these columns are only present when $matrix is present.
#           colidx      index of corresponding column in output matrix - this is the index of the simplified generators
#           cmax        coordinate with largest absolute value, used to compute the following
#           mixed       logical, which is TRUE  iff the group has "mixed directions"
#           major       the longer vector in the zonoseg spanned by the generators in the group, always non-zero
#           minor       the shorter vector in the zonoseg; this is non-zero iff the group has "mixed directions"
#                   The row index in multiplesupp  matches the list index of multiple.  We need a name for this index.
#   crossprods      3 x n(n-1)/2 matrix of all possible crossproducts after unitizing. Only present for simple rank 3.
#                   If i<j and the i'th and j'th generators are g_i and g_j,
#                   then the column of crossprods[] corresponding to i and j is filled with g_i X g_j unitized.
#   crossprodidx    integer LUT from hyperplane index to the column index of crossprods. Only present for simple rank 3.
#                   Used in getnormal.matroid()
#   hyperplaneidx   integer LUT from crossprods column index to hyperplane index. Only present for simple rank 3.
#                   this can be used to take a pair of points and find the unique hyperplane that contains them.
#                   Used in getmetrics.zonohedron()
#   collapsetosimple    integer map from original matrix columns to the simplified one - only present if the matroid is not simple
#   simplified      a matroid = the simplification of the original matroid - only present if the matroid is not simple already


#   x       a numeric matrix with 1,2, or 3 rows
#   e0      used when nrow(x) >= 1
#   e1      used when nrow(x) >= 2
#   e2      used when nrow(x) == 3
#
#   ground  integer vector, strictly increasing positive integers with length(ground) == ncol(x)
#

matroid.matrix <- function( x, e0=0, e1=1.e-6, e2=1.e-10, ground=NULL, ... )
    {
    #cat( "matroid.matrix\n" )
    #cat( "e0=", e0, "e1=", e1, "e2=", e2,  "ground=", ground, '\n' )

    time0   = gettime()

    ok  = is.numeric(x)  &&  (nrow(x)  %in%  1:3)  &&  (nrow(x) <= ncol(x))
    if( ! ok )
        {
        log_level( ERROR, "matrix x is invalid." )
        return(NULL)
        }

    if( is.integer(x) ) storage.mode(x) = 'double'

    elist   = list(e0,e1,e2)
    ok  = all( sapply(elist,is.numeric) )  &&  all( sapply(elist,length) == 1L )
    if( ! ok )
        {
        log_level( ERROR, "One of e0,e1,e2 is invalid." )
        return(NULL)
        }

    #   process the ground set
    if( is.null(ground) )
        {
        ground  = intfromchar( colnames(x) )

        if( is.null(ground) )
            ground  =  1:ncol(x)
        }

    if( length(ground) != ncol(x) )
        {
        log_level( ERROR, "ground is invalid, because the length=%d is incorrect. It must be %d.",
                                length(ground), ncol(x) )
        return(NULL)
        }

    if( ! all( 0 < diff(ground) ) )
        {
        log_level( ERROR,  "ground is invalid, because it is not in strictly increasing order." )
        return(NULL)
        }

    #   look for small columns in x;  these are the loops
    loopmask    = apply( x, MARGIN=2, function(vec) { max(abs(vec)) } ) <= e0
    loopraw     = which( loopmask )

    #   check for too many loops
    ok  = length(loopraw) <= ncol(x) - nrow(x)
    if( ! ok )
        {
        log_level( ERROR, "The matrix is rank-deficient.  |loops| = %d > %d = |columns| - |rows|.",
                            length(loopraw), ncol(x) - nrow(x) )
        return(NULL)
        }

    #   cleanup the colnames of x
    colnames(x) = as.character( ground )


    #   convert from matrix column indexes to ground indexes
    loop   = ground[ loopraw ]

    #   extract the nonloops, and the matrix of nonloops
    gnd.noloops = ground[ ! loopmask ]
    x.noloops   = x[ , ! loopmask, drop=FALSE ]


    #print( x.noloops )

    #   this group of statements works even when rank is 1

    xunit   = normalizeColumns( x.noloops ) #; print( xunit )

    grp = findColumnGroups( xunit, e1, oriented=FALSE ) #; print(grp)

    #   nonloop is a list of vectors, and not a vector
    #   the nonloop may contain many vectors that are singletons
    nonloop = setlistfromvec( grp, gnd.noloops )


    #   condenseMatrix() returns a list with the simplified matrix, and the data frame multiplesupp
    condata = condenseMatrix( x, ground, nonloop )

    x.simple    = condata$matrix        # this matrix has length(nonloop) columns.

    gnd.simple  = sapply( nonloop, function(v) { v[1] } )   # extract the first point in each group


    #   make integer vector that maps from the columns of x to the columns of x.simple
    #   if a column of x is 0, a loop, it maps to NA_integer_
    collapsetosimple    = rep( NA_integer_, length(loopmask) )

    collapsetosimple[ ! loopmask ]  = collapsetosimple( nonloop, gnd.noloops )


    lenvec      = lengths( nonloop )    #sapply( nonloop, length )
    multiple    = nonloop[ 2 <= lenvec ]

    lmdata      = list( loop=loop, multiple=multiple )

    issimple    = length(loop)==0  &&  length(multiple)==0

    if( nrow(x) == 1 )
        {
        #   almost trivial case rank=1
        #   the loops are the single hyperplane,
        #   and the non-loops are a single multiple group

        #out =  matroid.list( list(loop), ground=ground )
        out = matroid1( loop, gnd.noloops, ground, x, condata$multiplesupp )

        if( ! issimple )
            {
            out$collapsetosimple = collapsetosimple

            out$simplified  = matroid1( integer(0), gnd.simple[1], gnd.simple, x.simple )

            #   record the original loops and multiples as an attribute
            #   attr(out$simplified,"lmdata") = lmdata

            attr(out$simplified$hyperplane,"lmdata") = lmdata
            }

        return( out )
        }



    #   nrow(x) is 2 or 3

    if( length(nonloop) < nrow(x) )
        {
        log_level( ERROR, "The matrix is rank-deficient. rank = %d < %d = |rows|.",
                            length(nonloop), nrow(x) )
        return(NULL)
        }


    time1   = gettime()
    timelm  = time1 - time0


    if( nrow(x) == 2 )
        {
        #   case rank=2 is fairly short
        out = matroid2( loop, nonloop, ground, x, condata$multiplesupp )

        if( ! issimple )
            {
            out$collapsetosimple = collapsetosimple

            #   extract only the 1st point from each hyperplane
            hyperfirst  = lapply( nonloop, function(v) { v[1] } )
            # cat( "======\n" ) ;         print( nonloop ) ;             print( hyperfirst )

            #matrix_simp = condenseMatrix( x, out$ground, nonloop )

            out$simplified  = matroid2( integer(0), hyperfirst, gnd.simple, x.simple )

            #   record the original loops and multiples as an attribute
            #attr(out$simplified,"lmdata") = lmdata

            attr(out$simplified$hyperplane,"lmdata") = lmdata
            }

        return( out )
        }

    #########   nrow(x) is now 3, which is the difficult case       ######

    # cat( "Computing nontrivial hyperplanes:\n" )

    crossprods  = allcrossproducts( x.simple )
    crossprods  = normalizeMatrix( crossprods, 2L ) #; print( str(crossprods) )

    #   set bysize=TRUE, so that non-trivial hyperplanes are in descending order by size
    grp         = findColumnGroups( crossprods, e2, oriented=FALSE, bysize=TRUE )

    #   hypersnt is the number of non-trivial groups
    hypersnt    = max( grp )

    time2       = gettime()
    timenthp    = time2 - time1


    log_level( INFO, "In %dx%d matrix, found %d non-trivial hyperplanes.", nrow(x), ncol(x), hypersnt )

    if( 0 < hypersnt )
        {
        hyperplane  = vector( hypersnt, mode='list' )

        #   find the non-trivial hyperplanes, with more than 2 points
        pair    = allpairs( ncol(x.simple) )

        #   make table from nontrivial hyperplane index to column in crossprods
        crossprodidx_nontriv    = integer(hypersnt)

        #   make table from grouped column in crossprods to hyperplane index
        hyperplaneidx   = integer( ncol(crossprods) )

        warncount   = 0

        for( i in 1:hypersnt )
            {
            mask            = grp==i
            subpair         = pair[ mask, , drop=FALSE]     #; print( subpair )
            hyperraw        = fastunion( as.integer(subpair) )
            hyperplane[[i]] = gnd.simple[ hyperraw ]  # sort( unique( as.integer(subpair) ) ) ]

            m   = length( hyperplane[[i]] )

            if( nrow(subpair) != (m*(m-1L))/2L )
                {
                #   the vector grp[] is missing some pairs that should be in this hyperplane
                #  fill out the vector grp[] with *all* m*(m-1)/2 pairs in this hyperplane
                pairscomplete   = allpairs( hyperraw )
                pairidx         = .Call( C_pairindex, pairscomplete, ncol(x.simple) )

                if( length(pairidx) != (m*(m-1L))/2L )
                    {
                    log_level( FATAL, "Internal error.  length(pairidx)=%d  !=  %d.",
                                        length(pairidx), m*(m-1L)/2L )
                    return(NULL)
                    }

                grp[ pairidx ]  = i
                mask[ pairidx ] = TRUE      # retro-correct mask[] too !

                warncount   = warncount + 1

                if( warncount <= 10 )
                    {
                    log_level( WARN, "hyperplane %d came from %d pairs, but expected %d pairs.",
                                        i, nrow(subpair), m*(m-1L)/2L )
                    }
                }

            #   get normal from the first pair in the group
            whicheqi    = which(mask)
            crossprodidx_nontriv[i] = whicheqi[1]
            hyperplaneidx[whicheqi] = i

            #   normalnontriv[i, ] = cprods[ ,j]
            }

        if( 0 < warncount )
            log_level( WARN, "There were %d total warnings about hyperplanes and pairs.", warncount )

        #print( hyperplane )
        #print( x.simple )

        #time_start  = gettime()

        # pcount  = paircount1( hyperplane, gnd.simple )

        #time_elapsed    = gettime() - time_start
        #cat( "paircount1() time_elapsed=", time_elapsed, "sec\n" )

        #   find list of the trivials, with 2 points
        hypertriv   = trivialhypers2( hyperplane, gnd.simple )

        if( is.character(hypertriv) )
            {
            #   there was an error and hypertriv is the message
            hypertriv   = c( hypertriv, sprintf( "Try reducing argument e2=%g.", e2 ) )
            log_level( WARN, paste(hypertriv,collapse='\n') )
            #cat( paste(hypertriv,collapse='\n'), '\n' )
            return(NULL)
            }

        #print( hypertriv )

        idxtriv    = which(grp==0)

        crossprodidx_triv   = idxtriv
        crossprodidx        = c( crossprodidx_nontriv, crossprodidx_triv )

        hyperplaneidx[idxtriv]  = (1:length(idxtriv)) + hypersnt

        #   within each non-trivial hyperplane, snap all the cross products to agree exactly, up to sign
        # crossprods  = snapcrossprods( crossprods, hyperplane, crossprods[ ,crossprodidx_nontriv,drop=F], gnd.simple )
        # crossprods is modifed "in-place"
        .Call( C_snapcrossprods, crossprods, hyperplane, crossprods[ ,crossprodidx_nontriv,drop=F], gnd.simple )

        hyperplane  = c( hyperplane, hypertriv )    # concatenate the 2 lists

        if( length(hyperplane) != length(crossprodidx) )
            {
            log_level( FATAL, "Internal error.  nrow(normaltriv)=%d  !=  %d=length(hypertriv).",
                                length(hyperplane), length(crossprodidx) )
            return(NULL)
            }

        #   1st test
        ok  = 1 <= hyperplaneidx  &  hyperplaneidx <= length(hyperplane)
        if( ! all(ok) )
            {
            log_level( FATAL, "Internal Error.  %d values of hyperplaneidx are invalid.", sum(!ok) )
            return(NULL)
            }

        #   2nd test
        seqq    = 1:length(hyperplane)
        ok  = hyperplaneidx[ crossprodidx[seqq] ] == seqq
        if( ! all(ok) )
            {
            log_level( FATAL, "Internal Error.  %d values of hyperplaneidx o crossprodidx are invalid.", sum(!ok) )
            return(NULL)
            }

        #   The matrix pcount[,] has only 0s and 1s in the upper triangle
        #   We must add a trivial hyperplane for each 0.
        #   Count the number of 0s.
        #extra   = which( pcount==0L  &  row(pcount)<col(pcount), arr=T )
        #count   = nrow(extra)
        #cat( "Adding", count, "trivial hyperplanes...\n" )
        #hyperpair   = vector( count, mode='list' )
        #for( i in 1:count )
        #    hyperpair[[i]]  = gnd.simple[ extra[i, ] ]


        #cat( "Total is now", length(hyperplane), "hyperplanes.\n" )
        }
     else
        {
        #   *ALL* the hyperplanes are trivial, and have size 2
        #   the matroid is uniform
        #   the number of hyperplanes is the number of all the pairs = n*(n-1)/2  n is the number of columns in x.simple
        hyperplane      = matrix2list( allpairs(gnd.simple), 1L )

        crossprodidx    = 1:ncol(crossprods)

        hyperplaneidx   = 1:ncol(crossprods)

        #normal      = t( cprods )
        }

    if( FALSE )
        {
        #   transfer each row of matrix hypertriv to item in a list
        hyperpair  = vector( nrow(hypertriv), mode='list' )
        for( i in 1:nrow(hypertriv) )
            hyperpair[[i]] = hypertriv[i, ]
        }



    time3   = gettime()
    timethp = time3 - time2


    #   degeneracy check
    #   these 2 tests are equivalent, mathematically
    #   but perform them separately in case of unknown error

    if( length( hyperplane[[1]] ) == length(gnd.simple) )
        {
        log_level( ERROR, "The matrix is rank-deficient. rank = 2 < 3 = |rows|.  The first hyperplane is equal to the ground set (simplified)." )
        return(NULL)
        }

    if( length(hyperplane) == 1 )
        {
        log_level( ERROR, "The matrix is rank-deficient. rank = 2 < 3 = |rows|.  Only one hyperplane (simplified)." )
        return(NULL)
        }


    #   simplified matroid

    colnames(x.simple) = as.character( gnd.simple )

    if( ! issimple )
        #   record the original loops and multiples as an attribute of the hyperplane list
        attr(hyperplane,"lmdata") = lmdata

    simplified  = matroid3( hyperplane, integer(0), list(), gnd.simple, x.simple )


    #   add all the cross products.  n*(n-1)/2  n is the number of columns in x.simple
    simplified$crossprods   = crossprods

    #   add the lookup table from hyperplane index to the column of crossprods
    #   this assigns a normal vector to each hyperplane which effectively chooses
    #   one facet from the facet-pair of the zonohedron, namely the facet that
    #   has this normal as the outward pointing normal
    simplified$crossprodidx     = crossprodidx

    simplified$hyperplaneidx    = hyperplaneidx


    time4   = gettime()
    timesimplified  = time4 - time3


    #   unsimplified matroid
    if( issimple )
        {
        #   the matroid is *already* simple
        out = simplified
        }
    else
        {
        #   make output by unsimplification
        hyper_un = unsimplify( hyperplane, loop, multiple, gnd.simple )

        if( is.null(hyper_un) ) return(NULL)

        out = matroid3( hyper_un, loop, multiple, ground, x, condata$multiplesupp )

        out$collapsetosimple = collapsetosimple

        #   record the original loops and multiples as an attribute
        #   attr(simplified,"lmdata") = lmdata

        out$simplified  = simplified
        }

    time5   = gettime()
    timeunsimplified  = time5 - time4

    timetotal   = time5 - time0

    if( 0 )
    {
    timeother   =  timetotal -timelm-timenthp-timethp-timesimplified-timeunsimplified

    cat( "loops+multiples:          ", timelm * 1000, " msec\n" )
    cat( "non-trivial hyperplanes:  ", timenthp * 1000, " msec\n" )
    cat( "trivial hyperplanes:      ", timethp * 1000, " msec\n" )
    cat( "simplified matroid:       ", timesimplified * 1000, " msec\n" )
    cat( "unsimplified matroid:     ", timeunsimplified * 1000, " msec\n" )
    cat( "other:                    ", timeother * 1000, " msec\n" )
    cat( "Total:                    ", timetotal * 1000,  " msec\n" )
    }


    return(out)
    }


#   x       list of hyperplanes
#   ground  integer vector of ground set, only used when rank=1 and otherwise ignored

matroid.list <- function( x, ground=NULL, ... )
    {
    #cat( "matroid.list\n" )

    #   ensure all are in integer mode and sorted
    #x   = lapply( x, function(v) { sort.int(as.integer(v)) } )

    tmp = base::unlist( x )
    if( ! is.integer(tmp) )
        {
        log_level( ERROR, "x has non-integers." )
        return(NULL)
        }


    #   scrub distracting attributes, if any
    #   attr( x, "lmdata" )  = NULL

    n   = length(x)

    if( n <= 1 )
        {
        #   a rank=1 matroid, this is a special case
        #   hyperplanes have rank 0, so a hyperplane contains only loops
        #   and there can be at most 1 hyperplane
        if( is.null(ground) )
            {
            log_level( ERROR, "For a rank=1 matroid, ground cannot be NULL." )
            return(NULL)
            }

        ground  = as.integer(ground)

        if( ! all( 0 < diff(ground) ) )
            {
            log_level( ERROR, "ground is invalid, because it is not in strictly increasing order." )
            return(NULL)
            }

        if( n == 0 )
            loop    = integer(0)
        else
            loop    = x[[1]]

        ok  = is.integer(loop)  &&  subset1( loop, ground )
        if( ! ok )
            {
            log_level( ERROR, "The loops are invalid.  They must be integral and a subset of the ground set." )
            return(NULL)
            }

        nonloop = setdiff( ground, loop )
        if( length(nonloop) == 0 )
            {
            log_level( ERROR, "Every point in the ground set is a loop, which is invalid.  The rank==0." )
            return(NULL)
            }

        out = matroid1( loop, nonloop, ground )

        if( ! is_simple(out) )
            {
            out$simplified  = matroid1( integer(0), min(nonloop), min(nonloop) )

            #   record the original loops and multiples as an attribute
            attr(out$simplified$hyperplane,"lmdata") = list( loop=loop, multiple=list(nonloop) )
            }

        return( out )
        }


    if( ! is.null(ground) )
        {
        log_level( WARN, "For a rank>1 matroid, argument ground is ignored." )
        }

    ground  = fastunion( x )


    out = list()

    class( out )    = c( "matroid", class(out) )


    #   to test for ranks 2 and 3, first simplify
    hypersimple = simplify( x, ground )

    groundsimple    = fastunion( hypersimple )

    lenvec  = lengths( hypersimple )    # sapply( hypersimple, length ) is much slower
    minlen  = min( lenvec )

    if( minlen == 0 )
        {
        log_level( ERROR, "The simplified hyperplane list has an empty hyperplane, which is invalid." )
        return(NULL)
        }


    out$ground      = ground
    out$hyperplane  = x         #; print( str(out$hyperplane) )


    if( minlen == 1 )
        {
        #   this is rank 2
        #   check the *all* hyperplanes have size 1, and there are no duplicates
        ok  = all( lenvec == 1L )  &&  anyDuplicated( unlist(hypersimple) )==0
        if( ! ok )
            {
            log_level( ERROR, "A rank=2 matroid is detected, but the hyperplanes are invalid." )
            return(NULL)
            }
        out$rank        = 2L
        }
    else
        {
        #   check for rank=3
        #cat( "checking rank 3\n" ) ; flush.console()

        #cmat    = paircount1( hypersimple, groundsimple )

        # hypersimple must satisfy the paving axioms,
        # so hypertriv must be an empty list
        hypertriv   = trivialhypers2( hypersimple, groundsimple )

        if( is.character(hypertriv) )
            {
            #   there was an error and this is the message
            # cat( paste(hypertriv,collapse='\n'), '\n' )
            log_level( WARN, paste(hypertriv,collapse='\n') )
            return(NULL)
            }

        if( 0 < length(hypertriv) )
            {
            #   ERROR
            mess    = "The hyperplanes do not satisfy the paving matroid properties for rank=3."
            mess    = c( mess, "    There are %d point pairs that are in no hyperplane." )
            mess    = c( mess, "    One such point pair is %d,%d." )
            mess    = paste0( mess, sep='\n' )

            pair    = hypertriv[[1]]

            log_level( ERROR, mess, length(hypertriv), pair[1], pair[2] )

            return(NULL)
            }

        out$rank    = 3L
        }

    lmdata  = attr( hypersimple, "lmdata" )

    if( is.null(lmdata) )
        {
        out$loop        = integer(0)
        out$multiple    = list()
        }
    else
        {
        out$loop        = lmdata$loop
        out$multiple    = lmdata$multiple
        }

    out$multiplesupp    = data.frame( contiguous=is_contiguousgroup( out$loop, out$multiple, out$ground ) )

    if( ! is_simple(out) )
        {
        if( out$rank == 2L )
            {
            out$simplified  = matroid2( integer(0), hypersimple, groundsimple )
            }
        if( out$rank == 3L )
            {
            out$simplified  = matroid3( hypersimple, integer(0), list(), groundsimple )
            }

        #   record the original loops and multiples as an attribute
        attr(out$simplified$hyperplane,"lmdata")    = lmdata   # list( loop=out$loop, multiple=out$multiple )
        }


    #attr( hypersimple, "lmdata" )  = NULL        #   remove distraction


    return( out )
    }




#   loop    integer vector of loops, can be empty
#   nonloop integer vector of nonloops, cannot be empty
#   ground  ground set in ascending order
#
#   conditions, which are not checked:
#       loop and nonloop are disjoint
#
#   returns a matroid of rank 1
matroid1 <- function( loop, nonloop, ground, matrix=NULL, multiplesupp=NULL )
    {
    out = list()

    class( out )    = c( "matroid", class(out) )

    out$ground      = ground    # sort( c( loop, nonloop ) )
    out$hyperplane  = list(loop)
    out$rank        = 1L
    out$loop        = loop

    if( length(nonloop) == 1 )
        out$multiple    = list()
    else
        out$multiple    = list( nonloop )


    #   make supplementary data.frame with 1 column
    msupp   = data.frame( contiguous=is_contiguousgroup( loop, out$multiple, ground ) )

    if( ! is.null(matrix) )
        {
        if( ncol(matrix) != length(ground) )
            {
            log_level( ERROR, "%d != %d.", ncol(matrix), length(ground) )
            return(NULL)
            }

        out$matrix = matrix

        if( is.null(multiplesupp) ) multiplesupp    = emptymultiplesupp(1)

        if( nrow(multiplesupp) != length(out$multiple) )
            {
            log_level( ERROR, "multiplesupp  %d != %d.", nrow(multiplesupp), length(out$multiple) )
            return(NULL)
            }

        #   add matrix-related columns
        msupp = cbind( msupp, multiplesupp )
        }

    out$multiplesupp = msupp

    return( out )
    }

#   loop    integer vector of loops
#   nonloop list of integer vectors of nonloops, which define a partition of the nonloops
#   ground  integer vector for the ground set, in ascending order
#
#   conditions, which are not checked:
#       2 or more parts of the nontrivial partition, i.e. length(nonloop) >= 2
#       loops and nonloops are disjoint
#
#   returns a matroid of rank 2
matroid2 <- function( loop, nonloop, ground, matrix=NULL, multiplesupp=NULL )
    {
    out = list()

    class( out )    = c( "matroid", class(out) )

    out$ground  = ground    # fastunion( loop, nonloop )  # sort( c( loop, unique(unlist(nonloop)) ) )

    if( 0 < length(loop) )
        {
        #   add all loops to each set in the partition to get the hyperplanes
        myfun <- function( vec ) { sort.int( c(loop,vec) ) }
        out$hyperplane  = lapply( nonloop, myfun )
        }
    else
        out$hyperplane  = nonloop

    out$rank    = 2L

    out$loop    = loop

    sizevec         = lengths(nonloop)  #sapply( nonloop, length )
    out$multiple    = nonloop[ 2 <= sizevec ]

    #   make supplementary data.frame with 1 column
    msupp   = data.frame( contiguous=is_contiguousgroup( loop, out$multiple, ground ) )

    if( ! is.null(matrix) )
        {
        if( ncol(matrix) != length(ground) )
            {
            log_level( ERROR, "%d != %d.", ncol(matrix), length(ground) )
            return(NULL)
            }

        out$matrix = matrix

        if( is.null(multiplesupp) ) multiplesupp    = emptymultiplesupp(2)

        if( nrow(multiplesupp) != length(out$multiple) )
            {
            log_level( ERROR, "multiplesupp  %d != %d.", nrow(multiplesupp), length(out$multiple) )
            return(NULL)
            }

        #   add matrix-related columns
        msupp = cbind( msupp, multiplesupp )
        }

    out$multiplesupp = msupp

    return( out )
    }




#   matroid3()
#
#   unlike matroid1() and matroid2(), this one does very little processing.
#   It depends on the caller to do that.
#
#   hyperplane  the hyperplanes for a matroid of rank 3.
#               if the matroid is simple these sets form a 2-partition (of the ground set of the simple matroid)
#   loop        integer vector of loops
#   multiple    list of disjoint multiple groups
#   ground      integer vector for the ground set, in ascending order
#
#   conditions, which are not checked:
#       loops are disjoint from the remaining sets
#       the first point in each multiple group is in some hyperplane
#
#   returns a matroid of rank 3
#
#   note that if loop[] and multiple[] are empty, hyperplane can be used as is

matroid3 <- function( hyperplane, loop, multiple, ground=NULL, matrix=NULL, multiplesupp=NULL )
    {
    out = list()

    class( out )    = c( "matroid", class(out) )

    if( is.null(ground) )   ground = fastunion( hyperplane, loop, multiple )   # sort.int( unique( c( loop, unlist(multiple), unlist(hyperplane) ) ) )

    out$ground  = ground

    out$hyperplane  = hyperplane

    out$rank    = 3L

    out$loop    = loop

    out$multiple    = multiple

    #   make supplementary data.frame with 1 column
    msupp   = data.frame( contiguous=is_contiguousgroup( loop, multiple, ground ) )

    if( ! is.null(matrix) )
        {
        if( ncol(matrix) != length(ground) )
            {
            log_level( ERROR, "%d != %d.", ncol(matrix), length(ground) )
            return(NULL)
            }

        out$matrix = matrix

        if( is.null(multiplesupp) ) multiplesupp    = emptymultiplesupp(3)

        if( nrow(multiplesupp) != length(multiple) )
            {
            log_level( ERROR, "multiplesupp  %d != %d.", nrow(multiplesupp), length(multiple) )
            return(NULL)
            }

        msupp   = cbind( msupp, multiplesupp )  # add matrix-related columns
        }

    out$multiplesupp = msupp

    return( out )
    }


getsimplified.matroid <- function( x, ... )
    {
    if( is.null(x$simplified) )
        return(x)
    else
        return( x$simplified )
    }

unsimplify.matroid <- function( x, loop=NULL, multiple=NULL, ... )
    {
    return( matroid( unsimplify(x$hyperplane,loop,multiple,x$ground)  )  )
    }


is_uniform.matroid <- function( x )
    {
    if( x$rank == 1 )
        return( length(x$loop)==0 )

    if( ! is_simple(x) ) return(FALSE)

    return( all( lengths(x$hyperplane) == x$rank-1L ) )
    }

is_paving.matroid <- function( x )
    {
    if( x$rank == 3L )
        return( is_simple(x) )
    else if( x$rank == 2L )
        return( length(x$loop)==0 )
    else if( x$rank == 1L )
        return( TRUE )

    log_level( FATAL, "Internal error.  rank=%d", x$rank )

    return( NA )
    }

is_simple.matroid <- function( x )
    {
    return( length(x$loop)==0  &&  length(x$multiple)==0 )
    }


getground.matroid <- function( x )
    {
    return( x$ground )
    }

getloop.matroid <- function( x )
    {
    return( x$loop )
    }

getmultiple.matroid <- function( x )
    {
    return( x$multiple )
    }

gethyperplane.matroid <- function( x )
    {
    return( x$hyperplane )
    }

getmatrix.matroid <- function( x )
    {
    return( x$matrix )
    }

#   x           a simple matroid of rank 3, not checked
#   hyperidx    an integer m-vector of hyperplane indexes of x
#               if NULL then take this to be *ALL* the hyperplanes
#
#   returns a 3xm matrix with the "distinguished" normal of the hyperplanes in the rows

getnormal.matroid <- function( x, hyperidx, ... )
    {
    #if( ! is_simple(x)  ||  x$rank !=3 )
    #    {
    #    log.string( FATAL, "Internal error.  The matroid is invalid." )
    #    return(NULL)
    #    }

    #if( any( length(x$crossprodidx) < hyperidx ) )
    #    {
    #    cat( "getnormal(). hyperidx=", hyperidx, '\n' )
    #    return(NULL)
    #    }

    if( is.null(hyperidx) )
        out = x$crossprods[  , x$crossprodidx ]
    else
        out = x$crossprods[  , x$crossprodidx[hyperidx], drop=FALSE ]

    return( out )
    }




#   x   a matroid
#   idx a single integer, which is in the ground set
#
#   returns the index of the multiple[[]] group that contains idx,
#   and if there is none, then returns 0L
#   uses a brute force search, maybe optimize later
#
#   in case of error returns NULL

getmultipleindex.matroid <- function( x, idx )
    {
    ok  = is.integer(idx)  &&  length(idx)==1
    if( ! ok )
        {
        log_level( ERROR, "idx=%s is invalid.", as.character(idx) )
        return(NULL)
        }

    if( !( idx %in% x$ground ) )
        {
        log_level( ERROR, "idx=%d is invalid.", idx )
        return(NULL)
        }

    for( i in seq_len( length(x$multiple) ) )
        {
        if( idx %in% x$multiple[[i]] )  return(i)
        }

    return( 0L )
    }


#   x       a matroid
#   subs    a vector of integers representing subset of x$ground, or a list of such vectors
#
#   returns an integer vector with length = length of subs
#
#   if a set is NOT a subset of ground, the corresponding integer is NA_integer_

rank <- function( x, subs )
    {
    if( ! inherits( x, "matroid" ) )
        {
        log_level( ERROR, "x is not a matroid." )
        return(NULL)
        }

    if( ! is.list(subs) )
        subs    = list( as.integer(subs) )

    #   verify that all sets in subs are subsets of x$ground
    bad = ! .Call( C_issubset, subs, x$ground )
    if( any(bad) )
        {
        log_level( WARN, "%d of %d subsets are not a subset of ground.",
                        sum(bad), length(bad) )
        # cat( mess )
        # return( NULL )
        }

    names.saved = names(subs)

    if( ! is_simple(x) )
        subs    = .Call( C_simplifygeneral, subs, x$ground, x$loop, x$multiple )

    if( x$rank == 1L )
        {
        #  a special case
        out = lengths(subs)

        if( 1L < max(out) )
            {
            log_level( FATAL, "Internal error.  For rank 1 matroid, max(out) = %d > 1.", max(out) )
            return( NULL )
            }
        }
    else if( x$rank == 2L )
        {
        #  a special case
        out = pmin( lengths(subs), 2L )
        #   return( out )
        }
    else if( x$rank == 3L )
        {
        n   = length(subs)
        out = integer(n)

        #   find the non-trivial hyperplanes of 3 or points

        #   we do not need to simplfy to the non-trivials now,
        #   because C_anyissuperset uses the fact that the lengths of the hyperplanes are decreasing
        #   and can optimize it
        #hypersnt    = x$hyperplane[ 3 <= lengths(x$hyperplane) ]   this takes too long

        for( i in 1:n )
            {
            if( bad[i] )    next

            set = subs[[i]]

            if( length(set) <= 2L )
                {
                out[i]  = length(set)
                }
            else   # if( 0 < length(hypersnt) )
                {
                #   rank(set) is either 2 or 3, depending on whether set is a subset of a hyperplane
                test    =  .Call( C_anyissuperset, x$hyperplane, set, TRUE )   # .Call( C_issuperset, hypersnt, set )

                if( test )      # any(test) )
                    out[i] = 2L
                else
                    out[i] = 3L
                }
            #else
            #    #   length(set) >= 3 but the length of all hyperplanes is < 3
            #    out[i] = 3L
            }
        }
    else
        {
        log_level( ERROR, "rank(x)=%g != 3.", x$rank )
        return( NULL )
        }

    #   mark sets that are not subsets of ground with NA
    out[ bad ]  = NA_integer_

    names(out)  = names.saved

    return( out )
    }

is_independent <- function( x, subs )
    {
    if( ! is.list(subs) )
        subs    = list( as.integer(subs) )

    out = (rank(x,subs) == lengths(subs))

    names(out)  = names(subs)

    return( out )
    }

is_loop.matroid <- function( x, point )
    {
    #names.saved = names(point)
    #point   = as.integer(point)

    out = is.finite( match( point, x$loop ) )

    #   make points not in the ground set NA
    out[ is.na( match(point,x$ground) ) ] = NA

    names(out)  = names(point)      #names.saved

    return(out)
    }

#   lst     a list of integer vectors
charsummary <- function( lst )
    {
    if( length(lst) == 0 )    return( '{}' )

    out = ''

    if( length(lst) <= 8 )
        {
        for( vec in lst )
            out = c( out, sprintf( "{%s}", paste(vec,collapse=' ') ) )

        out = paste( out, collapse='  ' )

        if( nchar(out) <= 80 )
            return(out)
        }

    out = ''
    sizevec = unlist( lapply( lst, length ) )
    sizeunq = sort( unique(sizevec) )
    for( size in sizeunq )
        out    = c( out, sprintf( "  [%d-point: %d]", size, sum( sizevec==size ) ) )

    out = paste( out, collapse='' )

    return( out )
    }




print.matroid  <-  function( x, ... )
    {
    cat( "ground set:           ", length(x$ground), " points   {", paste(x$ground,collapse=' '), "}\n", sep='' )

    lmdata  = attr(x$hyperplane,"lmdata")       # attr(x,"lmdata")

    if( ! is.null(lmdata) )
        {
        #   this matroid is simple and derived from an "original" matroid
        for( vec in lmdata$multiple )
            {
            mess = paste( vec, collapse=' ' )
            mess = sprintf( "                    Point %d corresponds to the multiple group {%s} in the original matroid.\n",
                                    vec[1], mess )
            cat( mess )
            }
        }

    cat( "hyperplanes:          ", length(x$hyperplane), "   ", charsummary(x$hyperplane), '\n', sep='' )

    cat( "rank:                 ", x$rank, '\n', sep='' )

    cat( "loops:                ", length(x$loop), "   {", paste(x$loop,collapse=' '), "}", '\n', sep='' )

    cat( "multiple groups:      ", length(x$multiple), "   ", charsummary(x$multiple), '\n', sep='' )

    cat( "uniform:              ", is_uniform(x), '\n', sep='' )
    cat( "paving:               ", is_paving(x), '\n', sep='' )
    cat( "simple:               ", is_simple(x), '\n', sep='' )
    cat( "contiguous:           ", all(x$multiplesupp$contiguous), '\n', sep='' )


    if( ! is.null(x$matrix) )
        {
        mess    = sprintf( "This matroid is constructed from a %dx%d real matrix.\n", nrow(x$matrix), ncol(x$matrix) )
        cat( mess )
        if( ncol(x$matrix) <= 10 )
            print( x$matrix )
        }


    if( ! is.null(x$simplified) )
        {
        #   print the simplified matroid, and indent 4 spaces
        cat( '\n' )
        cat( "The summary of the simplified matroid is:\n" )
        mess = paste( "    ", capture.output(  print(x$simplified)  ), '\n', sep='' )
        cat( mess )
        }

    return( invisible(TRUE) )
    }




#   x       a list of subsets of a ground set, the hyperplanes of a matroid
#   ground  integer vector, if NULL computed from x$ground
#
#   finds all loops and groups of multiples
#   returns a new list of the same length with:
#       *) all loops removed
#       *) all multiples removed, except for the first point in each group

simplify.list <- function( x, ground=NULL, ... )
    {
    if( is.null(ground) )   ground = fastunion(x)

    lmdata  = loopsandmultiples( x, ground )       # ;  print( str(lmdata) )

    out = .Call( C_simplify, x, ground, lmdata$loop, lmdata$multiple )

    #   record the original loops and multiples as an attribute
    attr( out, "lmdata" ) = lmdata

    return( out )
    }


#   x       a list of subsets of a ground set, the hyperplanes of a matroid
#   ground  integer vector, if NULL computed from x$ground
#
#   finds all loops and groups of multiples
#   returns a new list of the same length with:
#       *) all loops removed
#       *) all multiples removed, except for the first point in each group

simplify_old.list <- function( x, ground=NULL  )
    {
    #   attr( x, "ground" ) = ground

    lmdata  = loopsandmultiples( x, ground )       # ;  print( str(lmdata) )

    remove  = lmdata$loop

    #   for each group of multiples, remove all except the first point
    for( idx in lmdata$multiple )
        remove  = c( remove, idx[-1] )

    if( length(remove) == 0 )   return(x)   # no change

    #   for each hyperplane, remove every point in the vector remove

    if( 1 )
        {
        for( i in seq_len(length(x)) )
            {
            hp  = x[[i]]
            idx = match( remove, hp, nomatch=0  )
            if( any( 0 < idx ) )    x[[i]]  = hp[ -idx ]
            }
        }
    else
        {
        #       this is actually SLOWER - TODO:  write a C version
        myfun   <- function( hp )
            {
            idx = match( remove, hp, nomatch=0  )

            if( any( 0 < idx ) )
                out = hp[ -idx ]
            else
                out = hp

            return(out)
            }

        x   = lapply( x, myfun )
        }

    #   record the original loops and multiples as an attribute
    attr( x,"lmdata") = lmdata

    return( x )
    }



#   x           a list of integer vectors, representing subsets of a ground set
#   loop        an integer vector, with all points disjoint from x
#   multiple    a list of multiples groups, each group must intersect the ground set in 1 point
#   ground      union of all points in x, in ascending order

unsimplify.list <- function( x, loop=NULL, multiple=NULL, ground=NULL, ... )
    {
    lmdata  =   attr(x,"lmdata")

    if( ! is.null(lmdata) )
        {
        if( is.null(loop) )     loop        = lmdata$loop
        if( is.null(multiple) ) multiple    = lmdata$multiple
        }
    else
        {
        if( is.null(loop) )     loop    = integer(0)
        if( is.null(multiple) ) multiple = list()
        }

    if( length(loop)==0  &&  length(multiple)==0 )
        {
        #   nothing to do
        attr(x,"lmdata")   = NULL       # ensure that "lmdata" is truly NULL
        return( x )
        }

    if( is.null(ground) )   ground  = fastunion(x)

    out = .Call( C_unsimplify, x, ground, loop, multiple )

    return( out )
    }







###########     helper functions   ############


#   hyperplane      a list of integer vector, defining subsets of a ground set
#                   Undocumented: it may also have an attribute "ground" = the ground set
#   ground          integer vector, if NULL computed from hyperplane$ground

#   returns a list with items:
#       loop        an integer vector of indexes of loops
#       multiple    a list of integer vectors, each of which is a group of multiples
#
#   loop        a point is a loop iff it appears in every hyperplane
#   multiple

loopsandmultiples   <- function( hyperplane, ground=NULL )
    {
    if( ! is.list(hyperplane) )
        {
        log_level( ERROR, "Argument hyperplane is not a list." )
        return(NULL)
        }

    out = list()

    m   = length(hyperplane)
    if( m == 0 )
        {
        #   no loops or multiples
        out$loop        = list()
        out$multiple    = list()
        return(out)
        }

    tmp = base::unlist( hyperplane )    #; cat( "tmp=", tmp, '\n' )
    if( ! is.integer(tmp) )
        {
        log_level( ERROR, "Argument hyperplane has non-integers." )
        return(NULL)
        }

    #   ground  = attr( hyperplane, "ground" )

    if( is.null(ground) )
        {
        #tmp     = unique(tmp)
        #ground  = sort( hyperplane )
        ground  = fastunion( hyperplane )
        }
    else
        {
        if( ! is.integer(ground) )
            {
            log_level( ERROR, "ground is non-integer." )
            return(NULL)
            }

        #   verify that ground is increasing
        if( ! all( 0 < diff(ground) ) )
            {
            log_level( ERROR, "ground is not strictly increasing." )
            return(NULL)
            }

        #   verify that tmp is a subset of ground
        if( ! subset1(tmp,ground) )
            {
            log_level( ERROR, "One of the hyperplanes is not a subset of ground." )
            return(NULL)
            }

        #   imax    = max( tmp, ground )
        }

    gmax    = ground[ length(ground) ]    # largest possible index

    maskg   = logical( gmax )
    maskg[ground]   = TRUE      # to be used below


    #   create the counters
    if( 1 )
        {
        res     = incidencedata( hyperplane, ground )
        if( is.null(res) )  return(NULL)

        #   incident[] an integer vector.  incident[i] = # of hyperplanes that contain point i
        #   hash[] a real vector depending on the incidence pattern of the point, hash[i] can be very large, so use real
        incident    = res$incident  #; print( incident )
        hash        = res$hash
        }
    else
        {
        #   first version too slow
        incident    = integer( gmax )   # incident[i] = # of hyperplanes that contain point i
        hash        = double( gmax )    # hash function of point i, hash[i] can be very large, so use double

        for( j in 1:m )
            {
            hp  = hyperplane[[ j ]]     # hp is the set of points in hyperplane j

            incident[hp]    = incident[hp]  + 1L
            hash[hp]        = hash[hp] + j^2            # a large signature, so collisions not likely
            }
        #cat( "hash=", hash, '\n' )
        }


    out$loop   = which( incident == m )

    out$multiple    = list()    # grow this list one at a time - slow and not good

    #   the first grouping only uses the hash function, and so it is only crude and approximate
    grp = grpDuplicated( matrix(hash,1,length(hash)), MARGIN=2 )
    if( all( grp==0 ) )
        {
        #   no multiples
        return( out )
        }
    #cat( "grp=", grp, '\n' )

    #   loops cannot be multiples, so zero them
    grp[ out$loop ] = 0L

    #   points not in the ground set cannot be multiples, so zero them
    grp[ ! maskg ]  = 0L

    #cat( "after zeroing loops and points outside ground set, grp=", grp, '\n' )

    grplist = grplistfromvec( grp )     #;   cat( "grplist:\n" ) ; print( grplist )

    pcount = 0  # of parallel groups
    for( idx in grplist )
        {
        #cat( "idx=", idx, '\n' )

        #   removing loops and points not in ground set may generate invalid groups
        if( length(idx) <= 1 )  next    # not a valid group

        #   make incidence matrix for columns taken from idx
        #   this is a subset of the full incidence matrix, and so saves a lot of memory
        if( 1 )
            mat = incidencematrix( hyperplane, ground, idx )
        else
            {
            mat = matrix( FALSE, m, length(idx) )
            for( j in 1:m ) { mat[j, ] =  idx  %in%  hyperplane[[ j ]] }
            }

        #  now compute the true multiple groups, using the full vector instead of a hash function
        grpsub      = grpDuplicated( mat, MARGIN=2 )    #;       cat( "grpsub=", grpsub, '\n' )
        grpsublist  = grplistfromvec( grpsub )
        for( idxsub in grpsublist )
            {
            #cat( "idxsub=", idxsub, '\n' )

            #   removing loops and points not in ground set may generate invalid groups
            if( length(idxsub) <= 1 )  next    # not a valid group

            pcount = pcount + 1

            #cat( "idxsub=", idxsub, "  pcount=", pcount, '\n' )
            out$multiple[[pcount]]  = idx[ idxsub ]
            }
        }

    return( out )
    }

#   x   matroid
#   W   invertible matrix - 2x2 or 3x3

lintransform.matroid <- function( x, W )
    {
    if( is.null(x$matrix) )
        {
        log_level( WARN, "matroid is not vectorial, so returning the matroid unchanged." )
        return(x)   # not generated from a matrix, so nothing can be done
        }

    if( length(W) == 1 )
        W = W * diag( x$rank )

    #   check that W is OK
    ok  = is.matrix(W)  &&  all( dim(W) == c(x$rank,x$rank) )
    if( ! ok )
        {
        log_level( ERROR, "matrix W is invalid." )
        return(NULL)
        }

    Winv    = try( solve(W), silent=TRUE )
    if( class(Winv)[1] == "try-error" )
        {
        log_level( ERROR, "matrix W is not invertible." )
        return(NULL)
        }


    #   just copy from x to out, and then make selective changes !
    out = x

    out$matrix  = W %*% x$matrix

    if( 0 < nrow(out$multiplesupp) )
        {
        out$multiplesupp$major  = x$multiplesupp$major %*% t(W)
        out$multiplesupp$minor  = x$multiplesupp$minor %*% t(W)
        }

    if( ! is.null(x$crossprods) )
        {
        #   transform crossprods, using Winv
        crossprods  = t(Winv)  %*%  x$crossprods

        out$crossprods  = normalizeMatrix( crossprods, 2L )
        }


    if( ! is.null(x$simplified) )
        {
        out$simplified$matrix  = W %*% x$simplified$matrix

        if( ! is.null(x$simplified$crossprods) )
            {
            #   transform crossprods, using Winv
            crossprods  = t(Winv)  %*%  x$simplified$crossprods

            out$simplified$crossprods  = normalizeMatrix( crossprods, 2L )
            }
        }

    return( out )
    }




#   loop        integer vector of loops
#   multiple    list of disjoint multiple groups
#   ground      integer vector for the ground set, in ascending order
#
#   returns a logical vector the same length as multiple[[]]

is_contiguousgroup <- function( loop, multiple, ground )
    {
    m   = length(multiple)

    out = logical(m)

    if( m == 0 ) return( out )

    # subtract loops from ground set
    if( 0 < length(loop) )
        {
        idx = match( loop, ground )
        gnd.noloops = ground[ -idx ]
        }
    else
        gnd.noloops = ground

    for( i in 1:m )
        out[i]  = is_contiguous( multiple[[i]], gnd.noloops )

    return( out )
    }


#   x           a simple matroid of rank 3.  Not checked it takes too long.
#   hypersub    integer m-vector of indexes of a hyperplane (in the simple matroid)
#   gen         a generator/point in the hyperplane (in the ground set of the simple matroid)
#   normal      mx3 matrix of normal vectors to the m hyperplanes given by hypersub
#
#   returns     3xm vector of m facet diameters

getdiametermatrix <- function( x, hypersub, pcube, gen, normal )
    {
    .Call( C_diametermatrix, x$hyperplane, hypersub, pcube, gen, x$ground, normal, x$matrix, x$crossprods )
    }

getbeltdata <- function( x, hypersub, pcube, gen, normal )
    {
    .Call( C_beltdata, x$hyperplane, hypersub, pcube, gen, x$ground, normal, x$matrix, x$crossprods )
    }


#   x       a matroid
#   idx     a vector of raw indexes in the simplified matroid, think of them as column indexes
#           they should be distinct (not verified)

#   returns a vector of raw indexes in the original matroid.
#   If an individual idx[k] is in a group, then that index expands to all indexes in the original.

liftrawindexes  <- function( x, idx )
    {
    if( ! inherits( x, "matroid" ) )
        {
        log_level( ERROR, "x is not a matroid." )
        return(NULL)
        }

    if( is_simple(x) )  return( idx )       # no change !

    if( is.null(x$collapsetosimple) )
        {
        log_level( ERROR, "x$collapsetosimple is NULL." )
        return(NULL)
        }

    mask    = logical( length(x$collapsetosimple) )

    mask[ x$collapsetosimple  %in%  idx ]   = TRUE

    return( which(mask) )
    }


#   returns vector with column indexes of all groups with mixed directions
getmixed.matroid <- function( x )
    {
    if( is.null(x$matrix) )         return(NULL)    # matroid did not come from a matrix

    if( nrow(x$multiplesupp)==0 )   return(integer(0))    #  no multiple groups

    mixed   = x$multiplesupp$mixed  # logical vector

    return( x$multiplesupp$colidx[mixed] )
    }



###########     UseMethod() functions   ############

matroid <- function( x, ... )
{
    UseMethod('matroid')
}

simplify <- function( x, ... )
{
    UseMethod('simplify')
}

getsimplified <- function( x, ... )
{
    UseMethod('getsimplified')
}

unsimplify <- function( x, ... )
{
    UseMethod('unsimplify')
}

is_simple <- function( x )
{
    UseMethod('is_simple')
}

is_uniform <- function( x )
{
    UseMethod('is_uniform')
}

is_paving <- function( x )
{
    UseMethod('is_paving')
}

getground  <- function( x )
{
    UseMethod('getground')
}

gethyperplane  <- function( x )
{
    UseMethod('gethyperplane')
}

#rank  <- function( x, subs )
#{
#    UseMethod('rank')
#}

#is_independent  <- function( x, subs )
#{
#    UseMethod('is_independent')
#}

is_loop  <- function( x, point )
{
    UseMethod('is_loop')
}

getmultipleindex <- function( x, idx )
{
    UseMethod('getmultipleindex')
}

getmixed <- function( x )
{
    UseMethod('getmixed')
}

getloop <- function( x )
{
    UseMethod('getloop')
}

getmultiple <- function( x )
{
    UseMethod('getmultiple')
}




##################      deadwood below      #######################################

#   x           a simple matroid of rank 3.  Not checked it takes too long.
#   hyperidx    index of a hyperplane (in the simple matroid)
#   gen         a generator/point in the hyperplane (in the ground set of the simple matroid)
#
#   parameters 2 and 3 define a zonogon facet, and a pair of antipodal edges of that zonogon
#
#   returns the vector from the midpoint of one edge to the midpoint of the antipodal edge

getdiameter <- function( x, hyperidx, gen )
    {
    #if( ! is_simple(x)  ||  x$rank !=3 )
    #    {
    #    log.string( FATAL, "Internal error.  The matroid is invalid." )
    #    return(NULL)
    #    }

    #   get all ground set points of the hyperplane
    hyper   = x$hyperplane[[hyperidx]]

    #   convert from ground set to raw index
    generatoridx    = match( hyper, x$ground )

    if( any( is.na(generatoridx) ) )
        {
        log_level( FATAL, "Internal error.  Hyperplane %g is not a subset of the ground set.", hyperidx )
        return(NULL)
        }

    k   = match( gen, hyper )
    if( is.na(k) )
        {
        log_level( FATAL, "Internal error.  Generator %g is not in hyperplane %g.", gen, hyperidx )
        return(NULL)
        }

    #   genidx  = generatoridx[k]

    #   reorder generatoridx so that genidx comes first
    generatoridx    = c( generatoridx[k], generatoridx[-k] )

    normal  = x$crossprods[  , x$crossprodidx[hyperidx] ]    #; print(normal)

    out = .Call( C_diametervector, generatoridx, normal, x$matrix, x$crossprods )

    return(out)
    }

