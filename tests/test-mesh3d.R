library( zonohedra )
# library( rgl )        not needed if we write rgl::as.mesh3d()


# library( logger )
# log_threshold( TRACE, namespace='zonohedra' )


options( width=160 )


testClassics <- function( )
    {
    n   = length(classics.genlist)

    zonolist    = vector( n, mode='list' )

    for( k in 1:n )
        {
        zonolist[[k]]   = zonohedron( classics.genlist[[k]] )
        }

    cat( '\n\n' )
    cat( "---------  testClassics()  --------------\n", sep='' )
    flush.console()

    timestart   = zonohedra:::gettime()

    for( k in 1:n )
        {
        cat( "-----------   k=", k, '   ---------------\n', sep='' )
        mesh = rgl::as.mesh3d( zonolist[[k]] )
        }

    timeelapsed = zonohedra:::gettime() - timestart

    cat( "    classics=", n, "  time=", timeelapsed, "sec\n", sep='' )
    flush.console()

    return(TRUE)
    }

testColorimetry <- function( )
    {
    zono    = zonohedron( colorimetry.genlist[[2]] )

    timestart   = zonohedra:::gettime()

    mesh = rgl::as.mesh3d( zono )

    timeelapsed = zonohedra:::gettime() - timestart

    cat( "    colorimetry 1nm:  time=", timeelapsed, "sec\n", sep='' )
    flush.console()

    return(TRUE)
    }


if( ! testClassics() )          stop( "testClassics() failed !" )

if( ! testColorimetry() )       stop( "testColorimetry() failed !" )


cat( "\nPassed all mesh3d() tests !\n", file=stderr() )

