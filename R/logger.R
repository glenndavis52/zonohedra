
g.options   <- list( stoponerror = TRUE                                  #   must be logical
                    )
    
#   put fn() between timestamp and the msg    
layout_mine <- structure(
    function(level, msg, namespace="zonohedra",
                                    .logcall = sys.call(), .topcall = sys.call(-1), .topenv = parent.frame())
        {
        # cat( "obj_addr()=", obj_addr( .topcall[[1L]] ), '\n' )
        # cat( "deparse1 =", deparse1( .topcall[[1L]] ), '\n' )
        
        fn  = deparse1( .topcall[[1L]] )
        
        paste0( attr(level, 'level'), ' [', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '] ', namespace, "::", fn, '(). ', msg )
        },
    generator = quote(layout_mine())
)

appender_mine <- structure(
    function(lines)
        {
        cat(lines, file = stderr(), sep = '\n' )
        
        #   test for STOP
        if( any( grepl("^(ERR|FATAL)",lines ) )  )
            {
            stop( "Stopping, because level is ERROR or FATAL.", call.=FALSE )
            }
        },
    generator = quote(appender_mine())
    )


######          deadwood below      #####################


if( FALSE )
{
formatter_mine <- structure(
    function(fmt, ..., .logcall = sys.call(), .topcall = sys.call(-1), .topenv = parent.frame())
        {
        out = eval(sprintf(fmt, ...), envir = .topenv)
        
        print( str( .topcall ) )       
        
        # paste0( sprintf("%s(). ", , eval(sprintf(fmt, ...), envir = .topenv), sep='' )
        out
        }, 
    generator = quote(formatter_mine())
    )
    
    

formatter_mine <-     function(fmt, ..., .logcall = sys.call(), .topcall = sys.call(-1), .topenv = parent.frame() )
    {
    out = eval(sprintf(fmt, ...), envir = .topenv) ; cat(out,'\n')
    
    #cat( deparse1( .topcall[[1L]] ), '\n' )    #  Inf recursion
    #cat( deparse1( .logcall[[1L]] ), '\n' )    #  Inf recursion too
    
    cat( "obj_addr()=", obj_addr( .topcall[[1L]] ), '\n' )

    sp =    sys.parents()
    
    for( k in sp )
        {
        cat( "k=", k, "  obj_addr()=", obj_addr( sys.call(k) ), '\n' )
        }


    # paste0( sprintf("%s(). ", , eval(sprintf(fmt, ...), envir = .topenv), sep='' )
    out
    }
    
    
make_formatter_mine <-  function()
    {
    out <- function(fmt, ..., .logcall = sys.call(), .topcall = sys.call(-1), .topenv = parent.frame())
        {
        out = eval(sprintf(fmt, ...), envir = .topenv)   ;  cat(out,'\n')
        
        cat( deparse( .topcall[[1]] ), '\n' )       
        
        # paste0( sprintf("%s(). ", , eval(sprintf(fmt, ...), envir = .topenv), sep='' )
        out
        }
        
    return( out )
    }
    
    
    
    
}