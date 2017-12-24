#define CHECK_KN
#include "RNM.hpp"
#include <cassert>

using namespace std;  

//===================================================================
//  définition des 6 types de base de tableaux a 1,2 et 3 paramètres
//===================================================================

typedef double R;
typedef KN<R> Rn;
typedef KN_<R> Rn_;
typedef KNM<R> Rnm;
typedef KNM_<R> Rnm_;
typedef KNMK<R> Rnmk;
typedef KNMK_<R> Rnmk_;

// deux fonctions pour l'initialisation des tableaux

R f(long i){return i;}    
R g(long  i){return -i;}


int main()
{
  const int n=8;  
  cout << "******* Hello World, this is RNM use! ****" << endl << endl;

cout<<"================================================"<<endl;
cout<<"  déclaration et initialisation d'un vecteur    "<<endl;
cout<<"================================================"<<endl;

  Rn a(n,f);             // initialisation avec $a_i=f(i)$
  Rn b(n),c(n);          // initialisation par défaut à 0

  b = 5.*a;              
  cout << " a      = " <<  a << endl;
  cout << " b =5.*a= " <<  b << endl;

cout<<"================================================"<<endl;
cout<<"  les opérations vectorielles optimisées        "<<endl;
cout<<"================================================"<<endl;

  c  = a + b;             cout <<"c  = a + b      =" << c << endl;
  c  = 5. *b + a;
  c  = a + 5. *b;
  c  = a - 5. *b;
  c  = 10.*a - 5. *b;
  c  = 10.*a + 5. *b;

  c += a + b;
  c += 5. *b + a;
  c += a + 5. *b;
  c += a - 5. *b;
  c += 10.*a - 5. *b;
  c += 10.*a + 5. *b;

  c -= a + b;
  c -= 5. *b + a;
  c -= a + 5. *b;
  c -= a - 5. *b;
  c -= 10.*a - 5. *b;
  c -= 10.*a + 5. *b;   cout <<"c -= 10.*a + 5. *b=" << c << endl;

cout<<"================================================"<<endl;
cout<<"  extraction de sous-vecteurs                   "<<endl;
cout<<"================================================"<<endl;

  Rn     u(20,f), v(20,g);       //initialisation
  Rn_    u10(u,SubArray(10,5));  //  la partie [5,5+10[ du tableau u

  cout << "u                    =" << u   << endl;
  cout << "u10(u,SubArray(10,5))=" << u10 << endl<<endl;

  cout << " u(SubArray(10))   ="  << u(SubArray(10))   << endl;
  cout << " u(SubArray(10,5)) ="  << u(SubArray(10,5)) << endl;
  cout << " u(SubArray(8,5,2))="  << u(SubArray(8,5,2)) << endl<<endl;

  cout << " v                      = " << v << endl;
  v(SubArray(10,5)) += u10;
  cout << "v(SubArray(10,5)) += u10  " << v << endl;

cout<<"================================================"<<endl;
cout<<"  déclaration et initialisation d'une matrice   "<<endl;
cout<<"================================================"<<endl;

  Rnm A(n+2,n);                 //  une matrice $n+2\times n$

  for (int i=0;i<A.N();i++)     // ligne 
     for (int j=0;j<A.M();j++)  // colonne
       A(i,j) = 10*i+j;  

  cout << "A=" << A << endl;

cout<<"================================================"<<endl;
cout<<"  extraction de lignes, colonnes                "<<endl;
cout<<"================================================"<<endl;

  cout << "A1j=A( 1 ,'.') = " << A( 1 ,'.') << endl;   // la ligne  1
  cout << "Ai3=A('.', 3 ) = " << A('.', 3 ) << endl;   // la colonne 3
  cout << " A(5,'.')[1] = " << A(5,'.')[1]  << endl;
  cout << " A(5,1)      = " << A(5,1)       << endl;
  cout << " A('.',5)(1) = " << A('.',5)(1)  << endl<<endl;


  Rn_   Ai3(A('.', 3 ));         //  la colonne 3 de la matrice
  Rn    CopyAi3(A('.', 3 ));     //  une copie de la colonne 3
  assert( & A(0,3) == & Ai3(0)); //  vérification des adresses

  cout << "CopyAi3 = " << CopyAi3 << endl;
  CopyAi3[3]=100;

  cout << "CopyAi3 = " << CopyAi3 << endl;
  cout << "    Ai3 = " << CopyAi3 << endl << endl;

cout<<"================================================"<<endl;
cout<<"  extraction de sous-matrices                   "<<endl;
cout<<"================================================"<<endl;

  cout << "A=" << A << endl;

  Rnm   S(A(SubArray(3),SubArray(3))); // la sous matrice 3x3
  Rn_   Sii(S,SubArray(3,0,3+1));      // la diagonale de la matrice S(sans copie)

  cout << "S=  A(SubArray(3),SubArray(3))) = " << S <<endl;
  cout << "Sii     = " << Sii <<endl;

  cout << " A(SubArray(3,2),SubArray(2,4)) = " << endl;
  cout <<   A(SubArray(3,2),SubArray(2,4)) << endl<<endl;

  A(SubArray(3,2),SubArray(2,4)) = -1;
  A(SubArray(3,2),SubArray(2,0)) = -2;
  cout << "A=" << A << endl;

cout<<"================================================"<<endl;
cout<<"  opérations avec matrices                      "<<endl;
cout<<"================================================"<<endl;

    b = 1;            //tous les éléments de b sont =1

//  Rn Ab(n+2) = A*b; //erreur ! 

    Rn Ab(n+2);
    Ab = A*b;

    cout << " A         =" << A  << endl;
    cout << " b         =" <<  b << endl;
    cout << " Ab = A*b  =" << Ab << endl;

cout<<"================================================"<<endl;
cout<<"  tableaux à trois indices                      "<<endl;
cout<<"================================================"<<endl;

   Rnmk B(3,4,5);

   for (int i=0;i<B.N();i++)          // première dimension
         for (int j=0;j<B.M();j++)    // deuxième dimension
             for (int k=0;k<B.K();k++)// troisième dimension
                B(i,j,k) = 100*i+10*j+k;

  cout << " B      = " << B << endl<<endl;
  cout << " B(1  ,2  ,'.')= " << B(1  ,2  ,'.') << endl;
  cout << " B(1  ,'.',3  )= " << B(1  ,'.',3  ) << endl;
  cout << " B('.',2  ,3  )= " << B('.',2  ,3  ) << endl;
  cout << " B(1  ,'.','.')= " << B(1,'.'  ,'.') << endl;
  cout << " B('.',2  ,'.')= " << B('.',2  ,'.') << endl;
  cout << " B('.','.',3  )= " << B('.','.',3  ) << endl;

  cout << " B(1:2,1:3,0:3)   = "
       << B(SubArray(2,1),SubArray(3,1),SubArray(4,0)) << endl<<endl;

cout<<"==========copie de sous-tableaux ==========" <<endl;

  Rnmk Bsub(B(FromTo(1,2),FromTo(1,3),FromTo(0,3)));
  cout << Bsub << endl<<endl;

  B(SubArray(2,1),SubArray(3,1),SubArray(4,0))  = -1;
   cout << " B      = " << B << endl<<endl;

  B(SubArray(2,1),SubArray(3,1),SubArray(4,0)) += -1;
  cout << " B      = " << B << endl<<endl;

  return 0;

}
