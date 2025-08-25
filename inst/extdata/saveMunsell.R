

saveMunsell <- function( pathin= "../inst/extdata/Munsell-Finland-1600.txt", pathout="../inst/extdata/Munsell1600.txt" )
    {
    requireNamespace( "munsellinterpol" )
    requireNamespace( "spacesXYZ" )

    linevec     = readLines( pathin )

    idx = grep( "^ABB", linevec[1:100] )   # ; print(idx)

    if( length(idx)==0 )    return(FALSE)

    linevec = linevec[ idx[1]:length(linevec) ]

    linesexpected   = 1600

    n   = length(linevec)

    if( n != linesexpected )
        {
        cat( "linesdata =", n, "   linesexpected =", linesexpected, '\n' )
        return( FALSE )
        }

    #   hue lookup
    huenum  = c( A=2.5, B=5, C=7.5, D=10, E=1.25, F=3.75, G=6.25, H=8.75 )  # ;  print( huenum )

    huename = c( RR="R", YR="YR", YY="Y", GY="GY", GG="G", BG="BG", BB="B", PB="PB", PP="P", RP="RP" )  #;  print( huename )

    munsell = character( n )

    for( i in 1:n )
        {
        line    = linevec[i]

        if( grepl( "NEUT", line ) )
            {
            value   = as.double( substr(line,5,7) ) / 100
            munsell[i]  = sprintf( "N%g/", value )
            next
            }

        hue1    = huenum[ substr(line,1,1) ]
        hue2    = huename[  substr(line,2,3) ]

        val     = as.double( substr(line,4,5) ) / 10
        chroma  = as.double( substr(line,6,7) )

        munsell[i]  = paste( hue1, hue2, val, '/', chroma, sep='', collapse='' ) #; print( munsell[i] )
        }

    # print( munsell )

    #munsell     = c( "N10/", munsell )

    #   adapt to Illuminant E
    C   = spacesXYZ::XYZfromxyY( c( spacesXYZ::standardxy( "C.NBS" ), 1) )

    CtoD65 = spacesXYZ::CAT( C, 'E' )

    XYZ     = munsellinterpol::MunsellToXYZ( munsell )
    XYZ     = spacesXYZ::adaptXYZ( CtoD65, XYZ )   #; print( XYZ )

    header  = character(0)
    header  = c( header, "#   XYZ data converted from README.txt" )
    header  = c( header, "#  https://sites.uef.fi/spectral/munsell-colors-glossy-all-spectrofotometer-measured/"  )
    header  = c( header, '#  adapted to Illuminant E for simplicity' )    
    header  = c( header, '' )
    
    write( header, file=pathout )
    
    utils::write.table( XYZ, file=pathout, append=TRUE )

    return(TRUE)
    }



