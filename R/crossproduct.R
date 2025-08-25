
crossproduct  <-  function( v, w )
    {
    out = c( v[2]*w[3] - v[3]*w[2], -v[1]*w[3] + v[3]*w[1], v[1]*w[2] - v[2]*w[1] )
    
    return( out )
    }
        


#   allcrossproducts() computes the 3D cross product for every pair of columns of a 3xN matrix
#
#   A   a  3xN matrix
#    
#   returns a 3 x N(N-1)/2 matrix with the column order the same as allpairs()

allcrossproducts <- function( A )
    {
    ok  = is.double(A) && is.matrix(A) && nrow(A)==3
    if( ! ok )
        return(NULL)

    out = .Call( C_allcrossproducts, A )

    return( out )
    }        
    
    
#   vec     positive integer scalar or vector.  
#           A scalar is interpreted as a length, and equivalent to the vector being 1:n
#
#   returns an integer matrix with dimension n*(n-1)/2, with rows being all pairs in standard order

allpairs <- function( vec )
    {
    n   = length(vec)
    if( n == 1L )
        n   = as.integer(vec)

    if( n < 2L )    return(NULL)
    
    out = .Call( C_allpairs, n )
    
    if( 1L < length(vec) )
        {
        #   transfer from 1:n to vec[]
        #   this can take over 1 msec, when n >= 320
        dim.saved   = dim(out)
        out         = vec[ out ]
        dim(out)    = dim.saved
        }
        
    return(out)
    }
    
#   crossprods      3 x N(N-1) matrix with unitized crossproduct pairs in the columns
#   hyperplane      list with integer vectors, defining non-trivial hyperplanes from the ground set, length M
#   crossprodsref   3 x M matrix, with the "reference" crossproducts in the columns, one for each hyperplane
#   ground          N vector of increasing integers forming the ground set

snapcrossprods  <- function( crossprods, hyperplane, crossprodsref, ground )
    {
    n   = length(ground)
    if( ncol(crossprods) != n*(n-1)/2 ) return(NULL)

    m   = length(hyperplane)
    if( ncol(crossprodsref) != m ) return(NULL)
    
    idxfromground   = idxfromgroundfun( ground )
    
    out = crossprods
    
    for( k in 1:m )
        {
        #   get the "reference" crossproduct
        cpref   = crossprodsref[ , k ]
        
        imax    = which.max( abs(cpref) )
        
        signref = sign( cpref[imax] )
        
        idxraw  = idxfromground[ hyperplane[[k]] ]
 
        #   copy cpref to all the corresponding columns, while changing the sign when appropriate
        idxpair = allpairs( idxraw )        
        
        for( i in 1:nrow(idxpair) )
            {
            j   = PAIRINDEX( idxpair[i,1], idxpair[i,2], n )    
            out[ ,j] = sign( out[imax,j] ) * signref * cpref
            }
        }
    
    return( out )
    }
    
    
    
###############     deadwood below      ####################################

#   crossproducts2() computes the 3D cross product for every pair of columns of a 3xN matrix
#
#   A   a  3xN matrix
#
#   returns a data.frame with 2 columns:
#       idx         column indexes j1 and j2
#       crossprod   the cross product of columns j1 and j2

crossproducts2 <- function( A )
    {
    ok  = is.double(A) && is.matrix(A) && nrow(A)==3
    if( ! ok )
        {
        return(NULL)
        }
        
    p12 = base::crossprod( A[1, ,drop=F], A[2, ,drop=F] ) #; print( str(p12) )
    p13 = base::crossprod( A[1, ,drop=F], A[3, ,drop=F] )
    p23 = base::crossprod( A[2, ,drop=F], A[3, ,drop=F] )
    
    d12 = p12 - t(p12)
    d13 = p13 - t(p13)
    d23 = p23 - t(p23)
    
    n   = ncol(A)
    
    if( requireNamespace( 'arrangements', quietly=TRUE ) )
        idx = arrangements::combinations(n,2)       # faster
    else
        idx = t( utils::combn(n,2) )   # matrix of pairs.  slower
    
    return( t( cbind( d23[idx], -d13[idx], d12[idx] ) ) )
        
    
    out = data.frame( row.names=1:nrow(idx) )
    
    out$idx         = idx
    out$crossprod   = cbind( d23[idx], -d13[idx], d12[idx] )
    
    return(out)
    }
    
    