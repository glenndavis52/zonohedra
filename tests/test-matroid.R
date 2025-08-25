require( zonohedra )
options( width=144 )   

testSimplify <- function()
    {
    path = system.file( "extdata/xyz1931.5nm.txt", package='zonohedra' )
    
    df  = read.table( path, header=T, stringsAsFactors=F ) #; print( str(df) )
    
    ground  = df$wavelength
    xyz5    = rbind( df$x, df$y, df$z )
    colnames( xyz5 )    = ground
    
    M   = matroid( xyz5, ground=ground ) #;   print(M)
    
    hypersimple = simplify( M$hyperplane, M$ground )
    
    hyper   = unsimplify( hypersimple, M$loop, M$multiple )
    
    if( ! identical( M$hyperplane, hyper ) )    return(FALSE)
    
    
    path    = system.file( "extdata/ciexyz31_1.csv", package='zonohedra' )
    
    df  = read.table( path, header=T, sep=',', stringsAsFactors=F ) #; print( str(df) )
    
    ground  = df$Wavelength
    xyz1    = rbind( df$x, df$y, df$z )
    colnames( xyz1 )    = ground
    
    M   = matroid( xyz1, e2=1.e-10, ground=ground ) #;   print(M)
    
    hypersimple = simplify( M$hyperplane, M$ground )
    
    hyper   = unsimplify( hypersimple, M$loop, M$multiple )
    
    if( ! identical( M$hyperplane, hyper ) )    return(FALSE)
    
    return(TRUE)
    }
    
    
if( ! testSimplify() )  stop( "testSimplify() failed !" )
    
cat( "Passed all matroid tests !\n", file=stderr() )    