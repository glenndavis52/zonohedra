

g.microbenchmark    = FALSE     # logical value, whether the package microbenchmark is loaded.  It must be unlocked.


.onLoad <- function( libname, pkgname )
    {
    #   unlockBinding( "g.microbenchmark", asNamespace('colorSpec') )   # asNamespace(pkgname) here generates a NOTE !

    g.microbenchmark    <<- requireNamespace( 'microbenchmark', quietly=TRUE )  #;  cat( "g.microbenchmark=", g.microbenchmark, '\n' )

    if( requireNamespace( 'logger', quietly=FALSE ) )
        {
        #   log_formatter( formatter_mine )
        log_formatter( logger::formatter_sprintf, namespace="zonohedra" )   # force sprintf(), even if glue is installed
        log_layout( layout_mine, namespace="zonohedra" )                    # put fn() between timestamp and the msg    
        log_appender( appender_mine, namespace="zonohedra" )                # maybe stop on ERROR or FATAL
        log_threshold( WARN, namespace="zonohedra" )                        # default is INFO
        }
    }



.onAttach <- function( libname, pkgname )
    {
    info    = library( help='zonohedra' )        #eval(pkgname) ?
    info    = format( info )
    mask    = grepl( "^(Version|Built)", info )     #Title
    info    = gsub( "[ ]+", ' ', info[mask] )
    # mess    = sprintf( "Attaching %s", pkgname )
    mess    = paste( c( "Package: zonohedra", "Author: Glenn Davis", info ), collapse='.  ' )   #; cat(mess)
    packageStartupMessage( mess )

    #initOptions()
    }
