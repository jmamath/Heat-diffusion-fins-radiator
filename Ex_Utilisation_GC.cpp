//
//  main.cpp
//  sfem_1
//
//  Created by jmamath on 11/6/17.
//  Copyright © 2017 Platestorm. All rights reserved.
//
#include <fstream>
#include <cassert>
#include <algorithm>
using namespace std;

#include <time.h>
inline double CPUtime(){
#ifdef SYSTIMES
    struct tms buf ;
    if (times(&buf) !=-1)
        return ((double)buf.tms_utime+(double)buf.tms_stime)/(long) sysconf(_SC_CLK_TCK) ;
    else
#endif
        return ((double) clock())/CLOCKS_PER_SEC ;
}

#define KN_CHECK
#include "RNM.hpp"
#include "GC.hpp"

typedef double R ;
class MatriceLaplacien1D: VirtualMatrice<R> { public:
    typedef VirtualMatrice<R>::plusAx plusAx ;
    MatriceLaplacien1D() {} ;
    void addMatMul(const KN_<R> & x, KN_<R> & Ax) const ;
    plusAx operator*(const KN<R> & x) const {return plusAx(this,x) ;}
} ;
void MatriceLaplacien1D::addMatMul(const KN_<R> & x, KN_<R> & Ax) const {
    int n= x.N(),n_1=n-1;
    double h=1./(n_1), d = 2/h, d1 = -1/h;
    
    for (int i=1;i< n_1; i++)
        Ax[i] = (x[i-1] +x[i+1]) * d1 + x[i]*d ;
    Ax[0]=x[0] ;
    Ax[n_1]=x[n_1] ;
}


int main(int argc,char ** argv) {
    typedef KN<double> Rn; typedef KN_<double> Rn_; typedef KNM<double> Rnm; typedef KNM_<double> Rnm_; {
        int n=10;
        Rnm A(n,n),C(n,n),Id(n,n); A=-1 ;
        C=0 ;
        Id=0 ;
        Rn_ Aii(A,SubArray(n,0,n+1)); Rn_ Cii(C,SubArray(n,0,n+1)); Rn_ Idii(Id,SubArray(n,0,n+1));
        for (int i=0;i<n;i++)
            Cii[i]= 1/(Aii[i]=n+i*i*i);
        Idii=1 ;
        cout << A;
        Rn x(n),b(n),s(n);
        for (int i=0;i<n;i++) b[i]=i;
        cout << "GradienConjugue preconditionne par la diagonale " << endl;
        
        GradienConjugue(A,C,b,x,n,1e-10);
        s = A*x;
        cout << " solution : A*x= " << s << endl;
        cout << "GradienConjugue preconditionnee par la identity " << endl;
        
        GradienConjugue(A,MatriceIdentite<R>(),b,x,n,1e-6);
        s = A*x;
        cout << s << endl;
        }
    {
        cout << "GradienConjugue laplacien 1D par la identity " << endl;
        int N=100;
        Rn b(N),x(N);
        R h= 1./(N-1);
        b= h;
        b[0]=0 ;
        b[N-1]=0 ;
        
        R t0 = CPUtime();
        GradienConjugue(MatriceLaplacien1D(),MatriceIdentite<R>() ,b,x,N,1e-5); cout << " Temps cpu = " << CPUtime() - t0<< "s" << endl;
        R err=0;
        for (int i=0;i<N;i++)
        {
            R xx=i*h;
            err= max(fabs(x[i]- (xx*(1-xx)/2)),err); }
        cout << "Fin err=" << err << endl; }
    return 0;
   
}




