


########    a few classic zonohedra from https://www.ics.uci.edu/~eppstein/junkyard/ukraine/ukraine.html    ######

#   returns a list of 3xN matrices, with class "genlist" prepended

makeClassics    <- function()
    {
    phi = (1 + sqrt(5))/2
    s2  = sqrt(2)
    s3  = sqrt(3)

    out = list()

    #   cube   3 generators
    W   = diag(3)
    shortname   = "C"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "cube"
    out[[shortname]] = W


    #   rhombic dodecahedron   4 generators
    W   = matrix( c(1,1,1, 1,-1,1, 1,1,-1, 1,-1,-1), nrow=3, byrow=FALSE )
    W   = reorderGenerators( W ) #; print( W )
    shortname   = "RD"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "rhombic dodecahedron"
    out[[shortname]] = W



    #   Bilinski dodecahedron   4 generators
    W   = matrix( c(1,phi,0,  phi,0,1,  0,1,phi,  -1,phi,0), nrow=3, byrow=FALSE )
    shortname   = "BD"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "Bilinski dodecahedron"
    out[[shortname]] = W


    #   rhombic icosahedron   5 generators
    W   = matrix( c(1,phi,0,  phi,0,1,  0,1,phi,  -1,phi,0,  phi,0,-1), nrow=3, byrow=FALSE )
    W   = reorderGenerators( W ) #; print( W )
    shortname   = "RI"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "rhombic icosahedron"
    out[[shortname]] = W


    #   rhombo-hexagonal dodecahedron   5 generators
    W   = matrix( c(0,s3,1,  s3,0,1,  0,-s3,1,  -s3,0,1,  0,0,2), nrow=3, byrow=FALSE )
    W   = reorderGenerators( W ) #; print( W )
    shortname   = "RHD"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "rhombo-hexagonal dodecahedron"
    out[[shortname]] = W


    #   rhombic triacontahedron   6 generators
    W   = matrix( c(1,phi,0,  phi,0,1,  0,1,phi,  -1,phi,0,  phi,0,-1,  0,-1,phi), nrow=3, byrow=FALSE )
    W   = reorderGenerators( W ) #; print( W )
    shortname   = "RT"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "rhombic triacontahedron"
    out[[shortname]] = W

    #   truncated octahedron   6 generators
    W   = matrix( c(1,1,0, 1,-1,0, 1,0,1, 1,0,-1, 0,1,1, 0,1,-1), nrow=3, byrow=FALSE )
    W   = reorderGenerators( W ) #; print( W )
    shortname   = "TO"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "truncated octahedron"
    out[[shortname]] = W


    #   truncated rhombic dodecahedron   7 generators
    W   = matrix( c(1,1,1,  1,1,-1,  1,-1,1,  1,-1,-1,  s3,0,0,  0,s3,0,  0,0,s3), nrow=3, byrow=FALSE )
    W   = reorderGenerators( W ) #; print( W )
    shortname   = "TRD"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "truncated rhombic dodecahedron"
    out[[shortname]] = W



    #   truncated cuboctahedron   9 generators
    W   = matrix( c(1,1,0, 1,-1,0, 1,0,1, 1,0,-1, 0,1,1, 0,1,-1, s2,0,0, 0,s2,0, 0,0,s2),  nrow=3, byrow=FALSE )
    W   = reorderGenerators( W ) #; print( W )
    shortname   = "TC"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "truncated cuboctahedron"
    out[[shortname]]    = W



    #   rhombic enneacontahedron    10 generators
    W   = c( 1,1,1, 1,-1,1, -1,1,1, -1,-1,1,
            0,phi,phi-1,    0,phi,1-phi,
            phi-1,0,phi,    phi-1,0,-phi,
            phi,phi-1,0,    phi,1-phi,0
            )

    W   = matrix( W, nrow=3, byrow=FALSE )   #; print(W)
    shortname   = "RE"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "rhombic enneacontahedron"
    out[[shortname]] = W


    #   rhombic hectotriadiohedron    12 generators
    r4  = 1 + s2

    W  = c( 1,  1, r4,
            1, -1, r4,
            1,  1, -r4,
            1, -1, -r4,
             r4, 1,  1,
             r4, 1, -1,
            -r4, 1,  1,
            -r4, 1, -1,
              1,  r4, 1,
             -1,  r4, 1,
              1, -r4, 1,
             -1, -r4, 1
            )

    W   = matrix( W, nrow=3, byrow=FALSE )   #; print(W)
    shortname   = "RH"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "rhombic hectotriadiohedron"
    out[[shortname]] = W




    #   truncated icosidodecahedron   15 generators
    W   = c( 1,phi,phi-1,   1,-phi,phi-1,
            1,-phi,1-phi,   1,phi,1-phi,
            phi,1-phi,      1,phi,1-phi,-1,
            phi,phi-1,-1,   phi,phi-1,1,
            phi-1,1,phi,    phi-1,-1,-phi,
            phi-1,1,-phi,   phi-1,-1,phi,
            2,0,0, 0,2,0, 0,0,2 )

    W   = matrix( W, nrow=3, byrow=FALSE )   #; print(W)
    W   = reorderGenerators( W ) #; print( W )
    shortname   = "TI"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "truncated icosidodecahedron"
    out[[shortname]] = W



    #   truncated small rhombicosidodecahedron       21 generators
    W   = c( 1,0,-phi,  1,0,phi,
            0,-phi,1,   0,phi,1,
            -phi,1,0,   phi,1,0,
            1,phi,phi-1,    1,-phi,phi-1,
            1,-phi,1-phi,   1,phi,1-phi,
            phi,1-phi,1,    phi,1-phi,-1,
            phi,phi-1,-1,   phi,phi-1,1,
            phi-1,1,phi,    phi-1,-1,-phi,
            phi-1,1,-phi,   phi-1,-1,phi,
            2,0,0, 0,2,0, 0,0,2 )

    W   = matrix( W, nrow=3, byrow=FALSE )
    shortname   = "TSR"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "truncated small rhombicosidodecahedron"
    out[[shortname]] = W

    class( out )    = c( "genlist", class(out) )

    return( out )
    }

makeColorimetry <- function()
    {
    out     = list()


    #   xyz at 5nm
    path    = "../inst/extdata/xyz1931.5nm.txt"
    W   = as.matrix(  read.table( path, sep=' ', header=T )  )
    rownames(W) = W[ ,1]
    W   = t( W[ ,2:4] )
    shortname   = "xyz1931.5nm"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "xyz at 5nm step (1931)"
    out[[shortname]] = W


    #   xyz at 1nm
    path    = "../inst/extdata/ciexyz31_1.csv"
    W   = as.matrix(  read.table( path, sep=',', header=T )  )
    rownames(W) = W[ ,1]
    W   = t( W[ ,2:4] )
    shortname   = "xyz1931.1nm"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "xyz at 1nm step"
    out[[shortname]] = W


    #   lms at 1nm
    path    = "../inst/extdata/lms2000.1nm.csv"
    df = read.table( path, header=T, sep=','  )   #; print( str(df) )
    W   = as.matrix( df[  , 2:ncol(df) ] )  #; print( str(data) )
    W   = t(W)
    W[ is.na(W) ] = 0
    colnames(W) = df[ ,1]
    shortname   = "lms2000.1nm"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "lms at 1nm step"
    out[[shortname]] = W

    #   xyz at 5nm, modified by Judd and Voss
    path    = "../inst/extdata/ciexyzjv.csv"
    df = read.table( path, header=T, sep=','  )   #; print( str(df) )
    W   = as.matrix( df[  , 2:ncol(df) ] )  #; print( str(data) )
    W   = t(W)
    W[ is.na(W) ] = 0
    colnames(W) = df[ ,1]
    shortname   = "ciexyzjv.5nm"
    attr( W, "shortname" ) = shortname
    attr( W, "fullname" ) = "xyz at 5nm step (1978)"
    out[[shortname]] = W

    class( out )    = c( "genlist", class(out) )

    return( out )
    }



#   matgen  3 x N matrix of generators
reorderGenerators <- function( matgen )
    {
    require( zonohedra )

    zono    = zonohedron( matgen )

    if( ! is_pointed(zono) )    return( matgen )   # no change

    normal  = supportingnormal0( zono )

    if( any( is.na(normal) ) )  return( NULL )  # should not happen

    rot3x3  = goodframe3x3( normal )

    genrot = crossprod( rot3x3, matgen )      # same as t(rot3x3) %*% matgen

    z   = genrot[3, ]

    #   these must be all positive
    if( ! all( 0 < z ) )    return( NULL )  # should not happen

    x   = genrot[1, ]
    y   = genrot[2, ]

    theta   = atan2( y, x )

    #   sort in cclockwise order
    perm    = order( theta )        #; print( perm )

    out     = matgen[ , perm]

    return( out )
    }


saveDatasets  <- function( .path="../data/zonohedra.rda" )
    {
    savevec = character(0)


    classics.genlist    = makeClassics()
    savevec = c( savevec, "classics.genlist" )

    colorimetry.genlist = makeColorimetry()
    savevec = c( savevec, "colorimetry.genlist" )


    ##  finally ready to save it
    save( list=savevec, file=.path, compress='xz' )   #     'xz'  'gzip'  FALSE

    return( invisible(TRUE) )
    }




