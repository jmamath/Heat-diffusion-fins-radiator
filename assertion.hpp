#ifndef ASSERTION_HPP_
#define ASSERTION_HPP_
//  to compile all assertion
//#define ASSERTION
// to remove all the assert  
//#define NDEBUG
#include <assert.h>
#include <stdlib.h>
inline void fassert(bool i) {if (!i)  { exit(1);}}
#define ASSERTION(i)  0
#ifdef ASSERTION 
#undef ASSERTION
#define ASSERTION(i) fassert(i)
#endif
// 4 small function to test if you are inside a segment  [a,b[,[a,b],]a,b],]a,b[
template<class T> inline bool  Inside_of(const T &  a,const T & i,const T & b)
 {return (a <= i) && ( i < b);}
template<class T> inline bool  Inside_oo(const T &  a,const T & i,const T & b)
 {return (a <= i) && ( i <= b);}
template<class T> inline bool  Inside_fo(const T &  a,const T & i,const T & b)
 {return (a > i) && ( i <= b);}
template<class T> inline bool  Inside_ff(const T &  a,const T & i,const T & b)
 {return (a > i)  && ( i < b);}
#endif
