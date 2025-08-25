

#include <cstddef>   /* for size_t */
#include <functional>
#include <map>
//  #include <unordered_map>        //  not any faster

//  flVecMap        is short for "fixed-length vector map"

//  flVecMap is a map with key being std::vector<T>.  All keys have the same length.
//              the value is int


template <typename T>
class flVecMap {
private:
    typedef std::map<  std::vector<T>, int >                    this_maptype ;
    typedef std::pair< typename this_maptype::iterator, bool >  this_pairtype;


    std::size_t     mN ;    // the length of every key in the map, enforced

public:

    nanotime_t      mTimeInsertion ;        //  nanoseconds, a 64-bit unsigned integer
    nanotime_t      mTimeVertices ;         //  nanoseconds, a 64-bit unsigned integer

    this_maptype    mMap ;  //  the map

    //  N       length of all the vectors
    //  count   expected size of the map when finished

    flVecMap( std::size_t N, std::size_t count=0 )    //  constructor
        {
        mN  = N ;

        mTimeInsertion  = 0 ;
        mTimeVertices   = 0 ;
        //  if( 0 < count )   mMap.reserve( 10*count );  // reserve() is defined for std::unordered_map, but not std::map
        }

    std::size_t     getN() const    { return mN ; }

    //  returns the 1-based index of x
    int getIndex( const T* x, std::size_t N, bool complement )
        {
        if( N != mN )   return(0);  //  invalid

        std::vector<T>  thevec( N ) ;

        //  initialize thevec from x.   Is there a faster way ?
        if( complement )
            for( unsigned int k=0 ; k<N ; k++ )  thevec[k]   = 1 - x[k];     // 0 -> 1  and  1 -> 0
        else
            for( unsigned int k=0 ; k<N ; k++ )  thevec[k]   = x[k];

        //  if thevec is *NOT* in the map, its index will be the current map size + 1
        int idx = int(mMap.size()) + 1;

        this_pairtype   res = mMap.insert( std::pair< std::vector<T>, int >(thevec,idx)  );

        if( res.second )
            {
            //  insertion successful, 1st appearance of this vector in the map
            //  Rprintf( "first insertion, with idx=%d\n", idx );
            }
        else
            {
            //  insertion failed, because this vector is already in the map
            //  ignore the value of idx assigned above
            //  res.first holds the (key,value) from mMap
            idx = res.first->second ;
            }

        return( idx );
        }
};


//  timing tests show that <unordered_map> takes the same time as plain <map>
#if 0
template <typename T>
struct std::hash< std::vector<T> > {
public:
    std::size_t operator() (const std::vector<T> & x) const
        {
        size_t ans = 0;

        for( int i=x.size()-1 ; 0<=i ; i-- )    ans += i*i * x[i] ;

        return ans;
        }
};
#endif
