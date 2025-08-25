#include <cstddef>   /* for size_t */
#include <unordered_set>
#include <functional>

#include "rcVec.h"

/* Needs global randbit, lshift, rshift defined */

namespace std {

template <typename T>
struct hash <rcVec<T> > {
public:
    size_t operator() (const rcVec<T> & x) const
    {   size_t ans = 0;
        int i ;
        for(i = x.len - 1; i >=0; i--)
            ans ^= ((hash<T>()(*(x.x + x.eltShift * i) )) ^ randbit.randbit) + (ans<<lshift) + (ans>>rshift) ;
        return ans;
    }
};

template <>
struct hash <CharSEXP > {
public:
    size_t operator() (const CharSEXP & x) const
    {
        return std::hash<char *>()(const_cast<char *>(CHAR(x.sexp)) );
    }
};


template <>
struct hash <Rcomplex > {
public:
    size_t operator() (const Rcomplex & x) const
    {
        size_t tmp;
        tmp = std::hash<double>()(x.r);
        return tmp ^ ( ( std::hash<double>()(x.i)  ^ randbit.randbit) + (tmp<<lshift) + (tmp>>rshift) );
    }
};


}

template <typename T>
class vecSetHash {  // a set with key being rcVec type
private:
    rcVec<T> aRC;
    typedef std::unordered_set<rcVec<T> > rcvSetType;
    std::pair<typename rcvSetType::iterator,bool> retPair; // not used
    rcvSetType rcvSet; // using operator< of rcVec<T>

public:
    //vecSetHash() {rcvSet.max_load_factor(0.125);}
    void duplicatedMat      (const T* x, const int* nrow, const int* ncol, int* const out, bool const byRow=true, bool const fromLast=false);
    void anyDuplicatedMat   (const T* x, const int* nrow, const int* ncol, int* const out, bool const byRow=true, bool const fromLast=false);
};

template <typename T>
void vecSetHash<T>::duplicatedMat (const T* x, const int* nrow, const int* ncol, int* const out, bool const byRow, bool const fromLast)
{
    /* put a logical vector of duplicated rows of numeric matrix x into out */
    if(byRow){
        aRC.eltShift = aRC.nVec = (int)(*nrow);
        aRC.vecShift = 1;
        aRC.len = (int)(*ncol);
    }else{
        aRC.eltShift = 1;
        aRC.vecShift = aRC.len = (int)(*nrow);
        aRC.nVec = (int)(*ncol);
    }

    rcvSet.clear();

    // unordered_set reserve: may or may not change bucket_count & rehashing
    rcvSet.reserve( aRC.nVec );

    // set insert: if not previously inserted, the .second of returned pair is true; otherwise false. the .first is an iterator for the (previously) inserted element, which is not used.
    if (fromLast) {
        aRC.x=const_cast<T*>(x) + ( byRow ? (*nrow)-1 : ((*ncol)-1)*(*nrow) );
        for(int i=aRC.nVec-1; i>=0; aRC.x -= aRC.vecShift)
            out[i--] = (int) !(rcvSet.insert( aRC ).second);
    }else {
        aRC.x=const_cast<T*>(x);
        for(int i=0; i<aRC.nVec; aRC.x += aRC.vecShift)
            out[i++] = (int) !(rcvSet.insert( aRC ).second);
    }

}


template <typename T>
void vecSetHash<T>::anyDuplicatedMat (const T* x, const int* nrow, const int* ncol, int* const out, bool const byRow, bool const fromLast)
{
    /* put a logical vector of duplicated rows of numeric matrix x into out */
    if(byRow){
        aRC.eltShift = aRC.nVec = (int)(*nrow);
        aRC.vecShift = 1;
        aRC.len = (int)(*ncol);
    }else{
        aRC.eltShift = 1;
        aRC.vecShift = aRC.len = (int)(*nrow);
        aRC.nVec = (int)(*ncol);
    }

    rcvSet.clear();

    // unordered_set reserve: may or may not change bucket_count & rehashing
    rcvSet.reserve( aRC.nVec );

    out[0] = 0; // result when no duplicates are found
    // set insert: if not previously inserted, the .second of returned pair is true; otherwise false. the .first is an iterator for the (previously) inserted element, which is not used.
    if (fromLast) {
        aRC.x=const_cast<T*>(x) + ( byRow ? (*nrow)-1 : ((*ncol)-1)*(*nrow) );
        for(int i=aRC.nVec-1; i>=0; aRC.x -= aRC.vecShift, --i)
            if( !(rcvSet.insert( aRC ).second) ) {
                out[0] = i + 1;
                break;
            }
    }else {
        aRC.x=const_cast<T*>(x);
        for(int i=0; i<aRC.nVec; aRC.x += aRC.vecShift, ++i)
            if( !(rcvSet.insert( aRC ).second) ){
                out[0] = i + 1;
                break;
            }
    }

}

#include <unordered_map>

template <typename T>
class vecMapHash {  // a map with key being rcVec type
private:
    rcVec<T> aRC;
    typedef std::unordered_map<rcVec<T>, int  > rcvMapType;
    std::pair<typename rcvMapType::iterator,bool> retPair;
    rcvMapType rcvMap; // using operator< of rcVec<T>

public:
    bool grpDuplicatedMat(const T* x, const int* nrow, const int* ncol, bool byRow, int* const out, int gcount[3] );
};


//  x       pointer to first entry in the matrix
//  *nrow   # of rows
//  *ncol   # of columns
//  byRow   if TRUE compare rows, otherwise columns
//  out     computed integer array, length is # of rows
//  gcount  group counts, total, singleton, and non-singleton
//
//  the function fills out[] so that:
//      out[i] == 0 means the i'th row is unique, and is in its own singleton group
//      out[i] == k means the i'th row is in non-singleton group k of equal rows;
//                  group k always has 2 more more rows in it.
//
//  return value = the number of groups, including the singletons

template <typename T>
bool
vecMapHash<T>::grpDuplicatedMat(const T* x, const int* nrow, const int* ncol, bool byRow, int* const out, int gcount[3] )
    {
    /* put a logical vector of duplicated rows of numeric matrix x into out */
    if( byRow )
        {
        aRC.eltShift = aRC.nVec = (int)(*nrow);
        aRC.vecShift = 1;
        aRC.len = (int)(*ncol);
        }
    else
        {
        aRC.eltShift = 1;
        aRC.vecShift = aRC.len = (int)(*nrow);
        aRC.nVec = (int)(*ncol);
        }

    rcvMap.clear();

    // unordered_map reserve: may or may not change bucket_count & rehashing
    rcvMap.reserve( aRC.nVec );

    int grpId = 0;      //  id of the current non-singleton group
    int dups = 0 ;      //  number of duplicates

    //  key     the aRC
    //  value   {the index of the row or column where the vector first appears} + 1   (the 1-based R index)
    //  map insert: if not previously inserted, the .second of returned pair is true; otherwise false. the .first is an iterator for the (previously) inserted element.

    aRC.x   = const_cast<T*>(x);
    for( int i=0; i<aRC.nVec ; aRC.x+=aRC.vecShift, i++ )
        {
        //  try to insert key=aRC and value=i+1
        retPair = rcvMap.insert( std::pair<rcVec<T>, int>(aRC,i+1)  );    //  note the i+1

        if( retPair.second )
            {
            //  insertion successful, 1st appearance of this vector
            //  the value of retPair.first->second is i+1
            out[i]  = 0 ;
            //Rprintf( "i=%d  first retPair.first->second=%d  out[i]=%d\n", i, retPair.first->second,  out[i] );
            }
        else
            {
            //  insertion failed, because this vector is already in the map

            //  find 0-based index of the first appearance of this vector
            //  we always have j < i
            int j = retPair.first->second - 1 ; //  note the -1

            if( out[j] == 0 )
                {
                //  2nd appearance of this vector, create a new group and record it in both places in out[]
                out[i]  = out[j] = ++grpId ;
                dups    += 2;
                //Rprintf( "i=%d  retPair.first->second=%d  grpId=%d  out[i]=%d\n", i, retPair.first->second, grpId, out[i] );
                }
            else
                {
                //  3rd appearance or higher, get grpId from proper location in out[] and copy it to out[i]
                out[i]  = out[j];
                dups    += 1;
                //Rprintf( "i=%d  third  retPair.first->second=%d  out[i]=%d\n", i, retPair.first->second,  out[i] );
                }
            }
        }

    //  aRC.nVec-dups   = # of singleton groups
    //  grpId           = # of non-singleton groups

    gcount[1]   = aRC.nVec-dups;    //  # of singleton groups
    gcount[2]   = grpId;            //  # of non-singleton groups
    gcount[0]   = gcount[1] + gcount[2];

    return( true );
    }
