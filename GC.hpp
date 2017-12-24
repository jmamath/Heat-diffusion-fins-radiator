// exemple de programmation du gradient conjugue precondiconnée
template<class R,class M,class P> 
int GradienConjugue(const M & A,const P & C,const KN_<R> &b,KN_<R> &x,
                    int nbitermax,double eps)
{
   int n=b.N();
   assert(n==x.N());
   KN<R> g(n), h(n), Ah(n), & Cg(Ah);  // on utilise Ah pour stocke Cg  
   g = A*x;  
   g -= b;// g = Ax-b
   Cg = C*g; // gradient preconditionne 
   h =-Cg; 
   R g2 = (Cg,g);   
   R reps2 = eps*eps*g2; // epsilon relatif 
   for (int iter=0;iter<=nbitermax;iter++)
     {      
       Ah = A*h;
       R ro =  - (g,h)/ (h,Ah); // ro optimal (produit scalaire usuel)
       x += ro *h;
       g += ro *Ah; // plus besoin de Ah, on utilise avec Cg optimisation
       Cg = C*g;
       R g2p=g2; 
       g2 = (Cg,g);
       if (g2 < reps2) { 
          cout << " convergence " ;
          cout.width(4);
          cout  << iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl; 
          return 1;// ok 
          }
       if (iter%100==0) {
          cout << "             " ;
          cout.width(4);
          cout <<  iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl; }

       R gamma = g2/g2p;       
       h *= gamma;
       h -= Cg;  //  h = -Cg * gamma* h       
     }
   cout << " Non convergence de la méthode du gradient conjugue " <<endl;
   return 0; 
}
template <class R> 
class MatriceIdentite: VirtualMatrice<R> { public:
 typedef typename VirtualMatrice<R>::plusAx plusAx;
 MatriceIdentite() {}; 
 void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const {   Ax+=x; } 
 plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);} 
}; 
