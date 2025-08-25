






#   A           a matrix [possibly with NAs or NaNs ?]
#   eps         difference tolerance, used to 'collapse' one row at a time
#   oriented    if FALSE, then 2 columns that differ only in sign are considered the same.
#   bysize      sort the groups in decreasing order by size; requires extra work

#   returns
#       groupidx    an integer vector with length(groupidx) = ncol(A)
#                   0 means this column is a trivial singleton group (most common)
#                   a positive integer g, means this column belongs to non-trivial group g
#                   a column with an NA is always in its own singleton cluster

findColumnGroups <- function( A, eps, oriented, bysize=FALSE )
    {
    if( ! oriented )
        {
        A   = conditionalAntipodal( A, eps/2, MARGIN=2 )
        
        if( is.null(A) )    return(NULL)
        }

    #   collapse each row, usually 3 of them
    Acollapsed  = array( NA_real_, dim=dim(A) )
    for( i in 1:nrow(A) )
        {
        Acollapsed[i, ]  = collapseGroups1D( A[i, ], eps=eps )
        }

    #   the next function is accelerated with C++
    out = grpDuplicated( Acollapsed, MARGIN=2 )
    
    if( bysize  &&  ! is.null(out) )
        out = relabelGrpIndexes( out )

    return(out)
    }
    

    
    
#   group   an non-negative integer vector, where positive integers indicate membership in a group
#
#   returns a vector so that the groups are in descending order by size
    
relabelGrpIndexes <- function( group )
    {
    n   = max( group )
    
    if( n <=1 )    return(group)       # all 0s or only 1 group, so no change
    
    member  = vector( n, mode='list' )
    
    for( i in 1:n ) member[[i]] = which( group==i )
    
    lenvec  = lengths(member)
    
    perm    = order( lenvec, decreasing=TRUE )
    
    out = group
    
    for( i in 1:n ) out[ member[[ perm[i] ]] ]  = i
        
    return( out )
    }



#   vec     vector of doubles
#   eps     small non-negative number

#   A _group_ is a maximal set of elements in vec
#   with adjacent differences all <= eps

#   the function modifies vec[] by replacing each value in a group
#   by the mean of that group.
#   Exception: if the group contains a single integer (possibly with repeats),
#   then each value in the group is replaced by that integer.
#
#   returns:  the modified vec

collapseGroups1D <- function( vec, eps )
    {
    ok  = is.numeric(vec)  &&  ! any(is.na(vec))
    if( ! ok )
        {
        log_level( ERROR, "Argument vec is invalid." )
        return(NULL)
        }
    
    ok  = is.numeric(eps)  &&  length(eps)==1  &&  0<=eps
    if( ! ok )
        {
        log_level( ERROR, "Argument eps is invalid." )        
        return(NULL)
        }
        
    if( length(vec)<=1  ||  eps==0 )  return(vec)     # no change

    #   sort vec in increasing order
    perm    = order(vec)
    out     = vec[perm]
    
    #   change vector out[] "in place"
    ok  = .Call( C_collapseGroups1D_R, out, eps )
    
    if( ! ok )  return(NULL)
    
    #   restore original order
    out[perm]   = out
    
    return(out)
    }
    
    
    
#   A       a numeric matrix    
#   eps     small positive number
#   MARGIN  1 (vectors are the rows) or 2 (vectors are the columns)
#   for each vector, search for the first number whose absolute value > eps
#   If that value is negative then apply antipodal, and otherwise the identity.
#   So in the returned matrix, in each vector the first "significant" value is positive
#    
conditionalAntipodal  <- function( A, eps, MARGIN )
    {
    ok  = is.double(A) && is.matrix(A)
    if( ! ok )
        {
        return(NULL)
        }
        
    ok  = is.double(eps) && length(eps)==1
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
    
    
    #   make a deep (non-shallow) copy of A, because C_conditionalAntipodal() modifies in-place
    out = duplicate(A)
    
    #   change matrix out[] "in place"
    ok  = .Call( C_conditionalAntipodal, out, eps, MARGIN )
    
    if( ! ok )  return(NULL)    
    
    return( out )
    }
    

    
    
duplicate <- function(x)
    {
    .Call(C_duplicateR, x)
    }
    
    
obj_addr <- function(x)
    {    
    .Call(C_obj_addr,x)
    }
    
    
    
############        deadwood below  ##################

#   too slow    
conditionalAntipodal1 <- function( A, eps )
    {
    myfun   <- function( vec )
        {
        idx = which( eps < abs(vec) )
        
        if( length(idx)==0  ||  0<vec[idx[1]] )
            return( vec )
        else
            return( -vec )
        }
    
    return( base::apply( A, MARGIN=2, myfun ) )
    }
    
#   too slow    
conditionalAntipodal2 <- function( A, eps )
    {
    A   = t(A)
    
    #   extract the first non-zero entry in each row
    first   = apply( A, 1, function(r) { r[ which(eps<abs(r))[1] ] } )
    #   first[ ! is.finite(first) ] = 0     # change NAs to 0
    A   = t( sign(first) * A     )
    
    return(A)
    }
