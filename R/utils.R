




###########     argument processing     ##############
#
#   A   a non-empty numeric NxM matrix, or something that can be converted to be one
#
#   Nmin    the minimum allowed number of rows
#
#   returns such a matrix, or NULL in case of error

prepareNxM  <-  function( A, M, Nmin=1 )
    {
    ok  = is.numeric(A) &&  M*Nmin<=length(A)  &&  (length(dim(A))<=2)  # &&  (0<M)

    ok  = ok  &&  ifelse( is.matrix(A), ncol(A)==M, ((length(A) %% M)==0)  )

    if( ! ok )
        {
        #print( "prepareNx3" )
        #print( sys.frames() )
        mess    = substr( paste0(as.character(A),collapse=','), 1, 10 )
        #arglist = list( ERROR, "A must be a non-empty numeric Nx3 matrix (with N>=%d). A='%s...'", mess )
        #do.call( log.string, arglist, envir=parent.frame(n=3) )
        #myfun   = log.string
        #environment(myfun) = parent.frame(3)

        Aname = deparse(substitute(A))

        #   notice hack with 2L to make log.string() print name of parent function
        #log.string( c(ERROR,2L), "Argument '%s' must be a non-empty numeric Nx%d matrix (with N>=%d). %s='%s...'",
        #                            Aname, M, Nmin, Aname, mess )

        #   in the next call, note the assignment to .topcall,
        #   which makes log_level() print the name of the calling function, and *not*  "prepareNxM()".
        #   Currently, this is the only place in the package where this is done.
        log_level( ERROR, "Argument '%s' must be a non-empty numeric Nx%d matrix (with N>=%d). %s='%s...'",
                                    Aname, M, Nmin, Aname, mess, .topcall=sys.call(-1L) )
        return(NULL)
        }

    if( ! is.matrix(A) )
        A = matrix( A, ncol=M, byrow=TRUE )

    return( A )
    }









#   A   an integer vector, representing a set of integers
#   B   an integer vector, representing a set of integers
#
#   returns  TRUE or FALSE, depending on whether A is a subset of B
#
#   subset1() is fastest, because there is no error checking

subset1 <- function( A, B )
    {
    all( is.finite( match(A,B) ) )
    }

subset2 <- function( A, B )
    {
    length( setdiff(A,B) ) == 0
    }


normalizeColumns <- function( mat )
    {
    normvec = sqrt( .colSums( mat^2, nrow(mat), ncol(mat) ) )

    return( t( t(mat) / normvec ) )
    }

#   A       a numeric matrix
#   MARGIN  1 (vectors are the rows) or 2 (vectors are the columns)
#   for each vector, divide by L2 norm to get a unit vector
#
normalizeMatrix  <- function( A, MARGIN )
    {
    ok  = is.double(A) && is.matrix(A)
    if( ! ok )
        {
        return(NULL)
        }

    MARGIN  = as.integer(MARGIN)
    ok  = length(MARGIN)==1  &&  MARGIN %in% 1L:2L
    if( ! ok )
        {
        return(NULL)
        }


    #   make a deep (non-shallow) copy of A, because C_normalizeMatrix() modifies in-place
    out = duplicate(A)

    #   change matrix out[] "in place"
    ok  = .Call( C_normalizeMatrix, out, MARGIN )

    if( ! ok )  return(NULL)

    return( out )
    }

#   ground  positive integer vector in strictly increasing order.  Not checked.
#
#   returns lookup table from these integers to the order in the sequence
#
#       out[idx]  is equivalent to match( idx, ground )  but faster

idxfromgroundfun <- function( ground )
    {
    #   make lookup table from ground to raw column index
    numgen  = length(ground)

    out = integer( ground[numgen] )
    out[ ground ] = 1:numgen

    return( out )
    }


shortlabel <- function( ivec )
    {
    if( length(ivec) <= 3 )
        out = paste(ivec,collapse='+')
    else
        {
        ran = range(ivec)
        out = sprintf( "%g+...+%g", ran[1], ran[2] )
        }

    return(out)
    }



#   x       a numeric matrix
#   e0      used when nrow(x) >= 1
#   e1      used when nrow(x) >= 2
#   ground  integer vector labeling the columns of x
#
#   the columns of x are taken as generators of a matroid
#
#   returns a matrix with no "loops" and "multiple groups"

simplify.matrix  <- function( x, e0=0, e1=1.e-6, ground=NULL, ... )
    {
    ok  = is.numeric(x)  &&  is.matrix(x)
    if( ! ok )
        {
        log_level( ERROR, "matrix x is invalid." )
        return(NULL)
        }

    if( is.integer(x) ) storage.mode(x) = 'double'

    elist   = list(e0,e1)
    ok  = all( sapply(elist,is.numeric) )  &&  all( sapply(elist,length) == 1L )
    if( ! ok )
        {
        log_level( ERROR, "One of e0,e1 is invalid." )
        return(NULL)
        }

    if( is.null(ground) )
        ground  =  1L:ncol(x)

    if( length(ground) != ncol(x) )
        {
        log_level( ERROR, "ground is invalid, because the length is incorrect." )
        return(NULL)
        }

    if( ! all( 0 < diff(ground) ) )
        {
        log_level( ERROR, "ground is invalid, because it is not in strictly increasing order." )
        return(NULL)
        }

    #   look for small columns in x;  these are the loops
    loopmask    = apply( x, MARGIN=2, function(vec) { max(abs(vec)) } ) <= e0

    loopraw     = which( loopmask )

    colnames(x) = as.character( ground )

    #   convert from matrix column indexes to ground indexes
    gnd.noloops = ground[ ! loopmask ]

    #   extract the matrix of nonloops
    x.noloops   = x[ , ! loopmask, drop=FALSE ]

    if( ncol(x.noloops) == 0 )  return(x.noloops)   # special case

    xunit   = normalizeColumns( x.noloops ) #; print( xunit )

    grp = findColumnGroups( xunit, e1, oriented=FALSE )   #; print(grp)

    if( all(grp==0) )   return(x.noloops)   # special case

    #   nonloop is a list of vectors, and not a vector
    nonloop = setlistfromvec( grp, gnd.noloops )

    condata = condenseMatrix( x, ground, nonloop )

    if( is.null(condata) )  return(NULL)

    return( condata$matrix )
    }


#   A           a numeric matrix
#   ground      integer vector labeling the columns of A
#   conspec     a list of integer vectors, all subsets of ground.
#               the corresponding columns of A are co-directional, with tolerance e1
#
#   for each vector in conspec take the dominant direction of the corresponding columns,
#   and copy that direction to output matrix.

#   returns a list with items:
#       matrix          with a column for each vector in conspec.  This is the simplified matrix.  It has length(conspec) columns.
#       multiplesupp    a data.frame with a row for each group of multiples in conspec, and these columns
#           colidx      index of corresponding column in output matrix - this is the index of the simplified generators
#           cmax        coordinate with largest absolute value, used to compute the following
#           mixed       logical, which is TRUE  iff the group has "mixed directions"
#           major       the longer vector in the zonoseg spanned by the generators in the group, always non-zero
#           minor       the shorter vector in the zonoseg; this is non-zero iff the group has "mixed directions"


condenseMatrix  <- function( A, ground, conspec )
    {
    ok  = is.double(A) && is.matrix(A)
    if( ! ok )
        {
        return(NULL)
        }

    ok  = length(ground) == ncol(A)
    if( ! ok )
        {
        return(NULL)
        }

    #   make inverse lookup vector
    colfromgnd  = integer( max(ground) )
    colfromgnd[ ground ]    = 1L:ncol(A)


    if( ! is.list(conspec) )
        {
        log_level( ERROR, "Argument conspec is invalid." )
        return(NULL)
        }

    if( ! subset1(unlist(conspec),ground))
        {
        return(NULL)
        }

    lenvec      = lengths( conspec )        #sapply( conspec, length )
    conspec1    = conspec

    # change all the non-singleton lists to any valid single value
    # these columns will be overwritten later
    multiplemask    = 2 <= lenvec

    conspec1[ multiplemask ]  = ground[1]

    out = list()

    out$matrix = A[ , colfromgnd[ unlist(conspec1) ], drop=F ]

    #   replace all the multiple columns with dominant direction
    colidx  = which( multiplemask )
    m       = length(colidx)

    cmax    = integer(m)
    major   = matrix( 0, nrow=m, ncol=nrow(A) )
    minor   = matrix( 0, nrow=m, ncol=nrow(A) )
    mixed   = logical(m)

    for( i in seq_len(m) )
        {
        k   = colidx[i]

        idx = colfromgnd[ conspec[[k]] ]

        res     = dominantDirection( A[ ,idx, drop=FALSE] )

        cmax[i] = res$cmax

        out$matrix[ ,k]    = res$dominant        # .rowSums( A[ ,idx], nrow(A), length(idx) )

        major[i, ]  = res$major
        minor[i, ]  = res$minor

        mixed[i]    = any( res$minor != 0 )
        }

    colnames(out$matrix)   = sapply( conspec, shortlabel )

    multiplesupp        = data.frame( row.names=colidx  )
    multiplesupp$colidx = colidx
    multiplesupp$cmax   = cmax
    multiplesupp$mixed  = mixed
    multiplesupp$major  = major
    multiplesupp$minor  = minor


    out$multiplesupp    = multiplesupp

    return( out )
    }


emptymultiplesupp   <-function( m )
    {
    out         = data.frame( row.names=character(0) )
    out$colidx  = integer(0)
    out$major   = matrix( 0, nrow=0, ncol=m )
    out$minor   = matrix( 0, nrow=0, ncol=m )
    out$mixed   = integer(0)

    return(out)
    }


#   nonloop     a list of integer vectors, with values in ground[]
#               the vectors should be a partition of ground[]
#   ground      an increasing vector of positive integers

collapsetosimple <- function( nonloop, ground )
    {
    #   make inverse lookup vector
    colfromgnd  = integer( max(ground) )
    colfromgnd[ ground ]    = 1:length(ground)

    out = rep( NA_integer_, length(ground) )

    for( i in 1:length(nonloop) )
        {
        idx = colfromgnd[ nonloop[[i]] ]

        out[idx] = i
        }

    return( out )
    }


#   A           a numeric matrix with rank 1
#
#   thus the columns of A are collinear and generate a zonoseg Z in space
#   there may be vectors in a single direction (0 is an endpoint of Z),
#   or both directions (0 is in the interior of Z)
#   In the latter case, we say the generators are "mixed".
#
#   returns a list with these items:
#       cmax        the 1-based index of the coordinate with maximum absolute value
#       dominant    the difference between the 2 endpoints of the zonoseg
#                   with the direction chosen so agree with the largest of sums in either direction
#                   NB:     out$dominant    = out$major - out$minor
#       major       the largest of the sums
#       minor       the smaller of the sums, most often it is 0, which means not "mixed"

dominantDirection   <- function( A )
    {
    ok  = is.double(A) && is.matrix(A)
    if( ! ok )
        {
        return(NULL)
        }

    #   find which row has the largest norm,
    #   since the rows are all collinear too, any norm will do so use L^inf
    cmax    = arrayInd( which.max(abs(A)), dim(A) )[1]
    rowmax  = A[ cmax, ]

    sumneg = rowSums( A[ , rowmax < 0, drop=FALSE ] )
    sumpos = rowSums( A[ , 0 < rowmax, drop=FALSE ] )

    out = list()

    out$cmax    = cmax

    delta   =  abs(sumpos[cmax]) - abs(sumneg[cmax])
    if( delta == 0 )
        #   this means that sumpos and sumneg are exact opposites
        #   in this case, choose the one with positive value at imax
        delta = sumpos[cmax] - sumneg[cmax]

    if( 0 < delta )
        {
        out$major = sumpos ;  out$minor = sumneg
        }
    else
        {
        out$major = sumneg ; out$minor = sumpos
        }

    out$dominant    = out$major - out$minor

    return(out)
    }


#   base        a non-zero n-vector
#   direction   a mxn matrix, whose rows are taken as n-vectors, AND are all multiples of base
#
#   returns     an m-vector of values, all are +1, -1, or 0, depending on the multiple
#
#   a halfspace test would work, but that would be a little slower

signvector <- function( base, direction )
    {
    if( length(base) != ncol(direction) )
        {
        log_level( FATAL, "Internal error. %d != %d.", length(base), ncol(direction) )
        return(NULL)
        }

    #   find which component of base has the largest norm,
    #   since the rows of direction are all collinear, any norm will do so use L^inf
    imax    =   which.max( abs(base) )

    out = sign( base[imax] ) * sign( direction[ , imax ] )

    return(out)
    }


#   grp     integer vector, as returned from grpDuplicated().  0s are the singletons.
#   ground  integer vector of the point indexes
#
#   returns a list  of integer vectors, as indexed by ground[],
#   each of which is a group of multiples with more than 1 point; compare with setlistfromvec().

grplistfromvec <- function( grp, ground=NULL )
    {
    if( is.null(ground) )
        ground = 1L:length(grp)
    else if( length(ground) != length(grp) )
        {
        log_level( ERROR, "length(ground)= %d is invalid!", length(ground) )
        return(NULL)
        }

    n   = max(grp)

    out = vector( n, mode='list' )
    for( i in seq_len(n) )
        out[[i]]    = ground[ which( grp == i ) ]

    return(out)
    }



#   grp     integer vector, as returned from grpDuplicated().  0s are the singletons.
#   ground  integer vector of the point indexes, with the same length as grp
#
#   returns a list of integer vectors, as indexed by ground[],
#   This list of sets includes both singletons and multiples; compare with grplistfromvec().
#   We always have length(out) <= length(grp)

setlistfromvec <- function( grp, ground=NULL )
    {
    if( is.null(ground) )
        ground = 1L:length(grp)
    else if( length(ground) != length(grp) )
        {
        log_level( FATAL, "length(ground)= %d is invalid!", length(ground) )
        return(NULL)
        }

    out = as.list( ground )

    keep    = ! logical( length(out) )

    n   = max(grp)
    for( i in seq_len(n) )
        {
        idx     = which( grp == i )
        out[[ idx[1] ]]    = ground[ idx ]

        #   remove all except the first point in the group
        keep[ idx[-1] ]   = FALSE
        }

    out = out[ keep ]

    return(out)
    }


#   xN  one of:
#       *) a list of integer vectors
#       *) an integer vector
#       *) NULL (ignored)
#
#   returns the union of all inputs, in strictly ascending order

fastunion <- function( x1, x2=NULL, x3=NULL )
    {
    return( .Call( C_fastunion, x1, x2, x3 ) )
    }



#   hyper      a list of integer vectors, each one of them defining a nontrivial hyperplane subset, and in increasing order.  Not checked
#   ground     an integer vector - the union of all the sets in hyper - and in increasing order.  Not checked

#   returns a list of 2-point hyperplanes that, combined with hyper, form a 2-partition of the ground set
#       if this is not possible, returns a character message with info on the problem

trivialhypers2 <- function( hyper, ground )
    {
    out = .Call( C_trivialhypers2, hyper, ground )

    if( ! is.null(out$cmax) )
        {
        #   ERROR.  make a good error message
        mess    = "The hyperplanes do not satisfy the paving matroid properties for rank=3."
        mess    = c( mess, "    the point pair %d,%d appears in %d hyperplanes." )
        mess    = paste0( mess, sep='\n' )
        mess = sprintf( mess, out$pmax[1], out$pmax[2], out$cmax )
        #   mess    = c( mess, "    Try reducing argument e2." )
        return( mess )
        }

    return( out )
    }

matrix2list <- function( x, MARGIN )
    {
    return( .Call( C_matrix2list, x, as.integer(MARGIN) ) )
    }

incidencedata <- function( hyper, ground=NULL )
    {
    if( is.null(ground) )   ground = fastunion(hyper)

    return( .Call( C_incidencedata, hyper, ground ) )
    }

incidencematrix <- function( hyper, ground, subset )
    {
    out =  .Call( C_incidencematrix, hyper, ground, subset )

    colnames(out)   = as.character(subset)

    return( out )
    }

findRunsTRUE <- function( mask, periodic=FALSE )
    {
    #   put sentinels on either end, to make things far simpler
    dif = diff( c(FALSE,mask,FALSE) )

    start   = which( dif ==  1 )
    stop    = which( dif == -1 )

    if( length(start) != length(stop) )
        {
        log_level( FATAL, "Internal error.  %d != %d", length(start), length(stop) )
        return(NULL)
        }

    stop    = stop - 1L

    if( periodic  &&  2<=length(start) )
        {
        m   = length(start)

        if( start[1]==1  &&  stop[m]==length(mask) )
            {
            #   merge first and last
            start[1]    = start[m]
            start   = start[ 1:(m-1) ]
            stop    = stop[ 1:(m-1) ]
            }
        }

    return( cbind( start=start, stop=stop ) )
    }


#   A   2xn matrix, with n vectors in the columns generating a convex cone
#       there are no 0s or multiples, so do not have to worry about generators in both directions
#
#   if the cone is salient, returns the indexes of the 2 columns, in CCW order
#   if the cone is not salient, returns integer(0)

cxconegenerators <- function( A )
    {
    n   = ncol(A)
    if( nrow(A)!=2  ||  n<=1 ) return(NULL)

    theta   = atan2( A[2, ], A[1, ] )
    perm    = order( theta )

    theta_sorted    = theta[perm]

    gapvec  = c( diff(theta_sorted), 2*pi - (theta_sorted[n] - theta_sorted[1]) )
    kmax    = which.max( gapvec )

    if( gapvec[kmax] <= pi )
        #   the vectors generate the entire plane, NOT salient
        return( integer(0) )

    #   there is a gap with angle > pi
    if( kmax < n )  k2 = kmax+1
    else            k2 = 1

    out = perm[ c(k2,kmax) ]         # ;match( c(k2,kmax), perm )

    return( out )
    }


#   u   numeric N-vector
#
#   returns Nx2 matrix with each row on the circle

tocircle <- function( u, tol=5.e-16 )
    {
    z   = exp( u * 1i )

    out = cbind( Re(z), Im(z) )

    idx = which( abs(out) < tol )

    if( 0 < length(idx) )   out[idx] = 0

    return(out)
    }


#   u, v    unit vectors of dimension n
#
#   returns nxn rotation matrix that takes u to v, and fixes all vectors ortho to u and v
#
#   Characteristic Classes, Milnor and Stasheff, p. 77

rotationshortest <- function( u, v )
    {
    n   = length(u)
    if( length(v) != n )    return(NULL)

    uv  = u + v
    if( all(uv == 0) )  return(NULL)

    out = diag(n)  -  (uv %o% uv)/(1 + sum(u*v))  +  2*(v %o% u)

    return(out)
    }

#   dir     non-zero 3-vector
#
#   returns a 3x3 orthogonal matrix where the 3rd column is parallel to dir
#   and the 1st and 2nd are orthogonal to dir
#   Thus frame3x3 rotates dir to the z-axis

goodframe3x3 <- function( dir )
    {
    out = base::svd( dir, nu=3 )$u

    #   the 1st column is a multiple of dir, but it might be a negative multiple
    if( sum( dir * out[ ,1] ) < 0 )
        # reverse sign of 1st column
        out[ ,1]  = -out[ ,1]

    #   move 1st column to the 3rd
    #   and ensure that determinant is positive, so it's a rotation
    if( 0 < det(out) )
        perm    = c(2,3,1)  # even
    else
        perm    = c(3,2,1)  # odd

    out = out[ , perm ]

    return( out )
    }



#   point   nx3 matrix of points on a great circle of S^2
#   normal  unit normal to the plane spanned by the great circle, this serves to orient the circle
#
#   returns a permutation of 1:n that puts the points in counter-clockwise order,
#   using the right-hand-rule with normal

orderoncircle <- function( point, normal )
    {
    #   find a 3x3 rotation matrix that rotates normal to the north pole
    pole    = c(0,0,1)

    if( all( abs(normal+pole) < 1.e-2 ) )
        {
        #   tweak normal away from -pole,
        #   this will change the projected circle to an ellipse but the order is still the same
        normal  = c(1,12,-12)/17     # from a Pythagorean quadruple
        }

    Q   = rotationshortest( normal, pole )  #; print( Q )

    #   compute xy, the z coordinate is very near 0 and ignored
    xy  = point %*% t( Q[1:2, ] )   #; print(xy) ; print( sqrt( rowSums(xy^2) ) )

    theta   = atan2( xy[ ,2], xy[ ,1] )

    perm    = order(theta)

    return( perm )
    }



#   normal  a non-zero 3-vector
#
#   returns a 3x2 matrix where the columns complete normal to an orthogonal basis
#           the 2 columns are unit vectors
#           if normal is added as a 3rd column, the 3x3 matrix preserves orientation
#           if normal is replaced by -normal, the columns of output matrix are swapped.

frame3x2fun <- function( normal, axischeck=FALSE )
    {
    ok  = is.numeric(normal)  &&  length(normal)==3  &&  0<sum(abs(normal))
    if( ! ok )
        {
        log_level( ERROR, "argument normal is invalid." )
        return(NULL)
        }

    out = base::svd( normal, nu=3 )$u[ , 2:3 ]     #; print( out )

    test    = crossproduct( out[ ,1], out[ ,2] )

    if( sum(test*normal) < 0 )
        {
        #   swap columns
        out = out[ , 2L:1L ] #; cat( 'frame3x2fun().  columns swapped !\n' )
        }

    if( axischeck &&  isaxis( -out[ ,1] ) )
        {
        out[ ,1]    = -out[ ,1]
        out         = out[ , 2L:1L ]    # swap
        }
    else if( axischeck && isaxis( -out[ ,2] ) )
        {
        out[ ,2]    = -out[ ,2]
        out         = out[ , 2L:1L ]    # swap
        }

    if( FALSE )
        {
        #   multiplication check
        test    = normal %*% out
        if( 1.e-14 < max(abs(test)) )
            {
            log_level( ERROR, "frame3x2fun() failed orthogonal test = %g!", max(abs(test))  )
            }

        test    = t(out) %*% out  -  diag(2)
        if( 1.e-14 < max(abs(test)) )
            {
            log_level( ERROR, "frame3x2fun() failed product test = %g!", max(abs(test))  )
            }
        }

    if( FALSE )
        {
        #   orientation check
        if( determinant( cbind(out,normal) )$sign <= 0 )
            {
            log_level( ERROR, "frame3x2fun() failed sign test.  det = %g!", det( cbind(out,normal) )  )
            }
        }

    return(out)
    }

isaxis <- function( vec )
    {
    n   = length(vec)

    return( sum(vec==0)==n-1  &&  sum(vec==1)==1 )
    }



#   vec     numeric vector
#   skip    integer coordinate to skip (ignore)

allequalskip <- function( vec, skip )
    {
    vec = vec[-skip]

    return( all( vec==vec[1] ) )
    }


gettime <- function()
    {
    return( microbenchmark::get_nanotime() * 1.e-9 )
    }


createtimer <- function()
    {
    now = microbenchmark::get_nanotime() * 1.e-9

    #   using a list
    return( list( created=now, now=now, elapsed=0, total=0 )  )

    #   using an array of doubles is actually slower !
    #out = c(now,now,0,0)
    #names(out) = c("created","now","elapsed","total")
    #return( out )
    }

updatetimer <- function( x, reset=TRUE )
    {
    now = microbenchmark::get_nanotime() * 1.e-9
    out = x

    out$elapsed = now - out$now
    out$total   = now - out$created
    if( reset ) out$now = now

    #   using integers
    #out[3L]  = now - out[2L]
    #out[4L]  = now - out[1L]
    #if( reset ) out[2L]  = now

    #   using names
    #out['elapsed']  = now - out['now']
    #out['total']  = now - out['created']
    #if( reset ) out['now']  = now

    return(out)
    }

#   vertex  2n x m matrix with vertices of a convex polygon in the rows
#
#   returns 4(n-1) x m matrix with tiling of the polygon into convex quadrilaterals
#           If m=3, it is ready to pass to rgl::quad3d()

makequads <- function( vertex )
    {
    n   = nrow(vertex) / 2L

    ok  = 2<=n  &&  as.integer(n)==n
    if( ! ok )
        {
        log_level( ERROR, "nrow(vertex)=%g is invalid.", nrow(vertex) )
        return(NULL)
        }

    mat = makequadindexes(n)    # matrix( c( 1:(n-1), 2:n, (2*n-1):(n+1), (2*n):(n+2) ), ncol=4 )

    idx = as.integer( t(mat) )      #; print(idx)

    out = vertex[ idx, ]

    return(out)
    }


#   n   # of generators of a zonogon, n>1
#       the vertices are indexed from 1 to 2n
#
#   returns (n-1) x 4 matrix with index-tiling of the zonogon into convex quadrilaterals
#               each row generates a quadrilateral

makequadindexes <- function( n )
    {
    cbind( 1:(n-1), 2:n, (2*n-1):(n+1), (2*n):(n+2) )
    }



#   p   a point in the n-cube, which we can think of as a transmittance spectrum
#
#   returns the min number of transitions between 0 and 1 necessary to achieve such a p, including interpolation
#   it always returns an even integer
#
#   this is the fast C version

transitioncount <- function( p )
    {
    .Call( C_transitioncount, p )
    }


#   p   a point in the n-cube, which we can think of as a transmittance spectrum
#
#   returns the min number of transitions between 0 and 1 necessary to achieve such a p, including interpolation
#   it always returns an even integer
#
#   this is the slow R version

transitioncount_old <- function( p )
    {
    n   = length(p)
    if( n == 0 )    return(0L)

    #   find all runs of points in the interior of [0,1], in a periodic way
    interior    = 0<p  &  p<1
    if( all(interior) )
        #   special case
        return( as.integer( 2 * floor( (n+1)/2 ) ) )

    mat     = findRunsTRUE( interior, periodic=TRUE )

    transitions = 0
    if( 0 < nrow(mat) )
        {
        inext   = c(2:n,1)
        iprev   = c(n,1:(n-1))

        for( i in 1:nrow(mat) )
            {
            start   = mat[i,1]
            stop    = mat[i,2]
            m       = stop - start + 1

            if( m < 0 ) m = m + n

            same    = p[ iprev[start] ] == p[ inext[stop] ]
            inc     = ifelse( same, 2*floor( (m+1)/2 ), 2*floor( m/2 ) )
            transitions = transitions + inc
            }
        }

    #   now remove all the interior coordinates
    ppure   = p[ ! interior ]

    #   add usual 0-1 transitions
    transitions = transitions  +  sum( diff(ppure) != 0 )

    #   add wrap-around transition, if present
    if( ppure[1] != ppure[ length(ppure) ] )    transitions = transitions + 1

    return( as.integer(transitions) )
    }


#   this is only valid if 1 <= i < j <= n
PAIRINDEX   <- function( i, j, n )
    {
    (i-1)*n - ((i)*(i+1))/2 + j
    }



#   subground   increasing integer M-vector, and subvector of ground, NOT verified
#   ground      increasing integer N-vector, with M <= N
#
#   each of these inputs has a matrix of all pairs, with standard ordering
#       for subground there are M*(M-1)/2 pairs
#       for ground there are N*(N-1)/2 pairs
#
#   each pair for subground also appears in as a pair for ground
#   The function returns an integer vector of length M*(M-1)/2
#   giving the row in ground of each pair from subground

translateallpairs <- function( subground, ground )
    {
    #   m   = length(subground)

    #subidxraw   = .Call( C_allpairs, m )
    #subidx      = subground[ subidxraw ]
    #dim(subidx) = dim(subidxraw)

    subidx  = allpairs( subground )

    #   subidx has ALL pairs from subground, in standard order
    #   each pair is also a pair from ground
    #   translate these to raw indexes in ground

    idxfromground   = idxfromgroundfun( ground )
    idxraw          = idxfromground[ subidx ]
    dim(idxraw)     = dim(subidx)

    n   = length(ground)

    #   idxraw has some pairs taken from 1:n
    #   find their position in the standard order of ALL pairs - N(N-1)/2 of them
    out     = .Call( C_pairindex, idxraw, n )

    return( out )
    }


#   .vec1 and .vec2     non-zero vectors of the same dimension
#
angleBetween  <-  function( .vec1, .vec2, unitized=FALSE, eps=5.e-14 )
    {
    q   = sum( .vec1*.vec2 )

    if( ! unitized )
        {
        len1    = sqrt( sum(.vec1^2) )
        len2    = sqrt( sum(.vec2^2) )     #;    print( denom )

        denom   = len1 * len2

        if( abs(denom) < eps )    return( NA_real_ )

        q   = q / denom  #; print(q)
        }

    if( abs(q) < 0.99 )
        {
        #   the usual case uses acos
        out = acos(q)
        }
    else
        {
        #   use asin instead
        if( ! unitized )
            {
            .vec1   = .vec1 / len1
            .vec2   = .vec2 / len2
            }

        if( q < 0 ) .vec2 = -.vec2

        d   = .vec1 - .vec2
        d   = sqrt( sum(d*d) )

        out = 2 * asin( d/2 )

        if( q < 0 ) out = pi - out
        }

    return(out)
    }



#   returns NULL unless ALL the values are valid integers

#   if some are invalid, prints a warning message

intfromchar <- function( charvec )
    {
    if( is.null(charvec)  ||  ! is.character(charvec) ) return(NULL)

    bad = ! grepl( "[ ]*-?[0-9.]+[ ]*", charvec )
    if( any(bad) )
        {
        log_level( WARN, "%d of %d values are invalid integers.", sum(bad), length(bad) )
        return(NULL)
        }

    out = as.integer( charvec )

    bad = (out != as.double(charvec))

    bad[ is.na(bad) ]   = TRUE

    if( any(bad) )
        {
        log_level( WARN, "%d of %d values are invalid integers.", sum(bad), length(bad) )
        return(NULL)
        }

    return( out )
    }


#################################       deadwood below      ##########################

#   A           a numeric matrix
#   ground      integer vector labeling the columns of A
#   conspec     one of these:
#                   *) integer vector defining a subset of ground
#                   *) a list of integer vectors, all subsets of ground
#   For the integer vector, just copy those columns to output.
#   For a list, for each vector take the dominant direction of the corresponding columns,
#   and copy that direction to output.

#   returns a matrix with length(conspec) columns

condenseMatrix_old  <- function( A, ground, conspec )
    {
    ok  = is.double(A) && is.matrix(A)
    if( ! ok )
        {
        return(NULL)
        }

    ok  = length(ground) == ncol(A)
    if( ! ok )
        {
        return(NULL)
        }

    #   make inverse lookup vector
    colfromgnd  = integer( max(ground) )
    colfromgnd[ ground ]    = 1L:ncol(A)

    if( is.integer(conspec) )
        {
        #   this is the easy case
        if( ! subset1(conspec,ground) )
            {
            return(NULL)
            }

        out = A[ , colfromgnd[conspec], drop=F ]
        }
    else if( is.list(conspec) )
        {
        if( ! subset1(unlist(conspec),ground))
            {
            return(NULL)
            }

        lenvec      = lengths( conspec )        #sapply( conspec, length )
        conspec1    = conspec

        # change all the non-singleton lists to any valid single value
        multiplemask    = 2 <= lenvec

        conspec1[ multiplemask ]  = ground[1]

        out = A[ , colfromgnd[ unlist(conspec1) ], drop=F ]

        #   replace all the multiple columns with sums
        multiple    = which( multiplemask )
        for( k in multiple )
            {
            idx = colfromgnd[ conspec[[k]] ]

            out[ ,k]    = dominantDirection( A[ ,idx, drop=FALSE] )$dominant  # .rowSums( A[ ,idx], nrow(A), length(idx) )
            }
        }
    else
        {
        log_level( ERROR, "conspec is invalid." )
        return(NULL)
        }

    colnames(out)   = sapply( conspec, shortlabel )

    return( out )
    }