// vérification d'allocation
#include "assertion.hpp" 
#include <cstdlib>
#include <math.h>
#include <fstream>     

// quelques fonctions utiles
template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}
template<class T> inline T Abs (const T &a){return a <0 ? -a : a;}
template<class T> inline double Norme (const T &a){return sqrt(a*a);}
template<class T> inline void Exchange (T& a,T& b) {T c=a;a=b;b=c;}
template<class T> inline T Max (const T &a,const T & b,const T & c){return Max(Max(a,b),c);}
template<class T> inline T Min (const T &a,const T & b,const T & c){return Min(Min(a,b),c);}

using namespace std;
// définition de R
typedef double R;

// La classe R2
class R2 {
  friend ostream& operator <<(ostream& f, const R2 & P )
  { f << P.x << ' ' << P.y   ; return f; }
  friend istream& operator >>(istream& f,  R2 & P)
  { f >>  P.x >>  P.y  ; return f; }
  
public:  
  R x,y;
  R2 () :x(0),y(0) {}
  R2 (R a,R b):x(a),y(b)  {}
  R2 (R2 a,R2 b):x(b.x-a.x),y(b.y-a.y)  {}
  R2   operator+(R2 P)const   {return R2(x+P.x,y+P.y);}
  R2   operator+=(R2 P)  {x += P.x;y += P.y;return *this;}
  R2   operator-(R2 P)const   {return R2(x-P.x,y-P.y);}
  R2   operator-=(R2 P) {x -= P.x;y -= P.y;return *this;}
  R2   operator-()const  {return R2(-x,-y);}
  R2   operator+()const  {return *this;}
  R   operator,(R2 P)const  {return  x*P.x+y*P.y;} // produit scalaire
  R   operator^(R2 P)const {return  x*P.y-y*P.x;} // produit mixte
  R2   operator*(R c)const {return R2(x*c,y*c);}
  R2   operator/(R c)const {return R2(x/c,y/c);}
  R2   perp() {return R2(-y,x);} // la perpendiculaire
};
R2 operator*(R c,R2 P) {return P*c;}


// La classe Label (références de sommets ou triangles) 
class Label {  
  friend ostream& operator <<(ostream& f,const Label & r  )
  { f <<  r.lab ; return f; }
  friend istream& operator >>(istream& f, Label & r  )
  { f >>  r.lab ; return f; }
public: 
  int lab;
  Label(int r=0):lab(r){}
  int onGamma() const { return lab;} 
};


// La classe Vertex (modélisation des sommets) 
class Vertex : public R2,public Label {
public:
  friend ostream& operator <<(ostream& f, const Vertex & v )
  { f << (R2) v << ' ' << (Label &) v   ; return f; }
  friend istream& operator >> (istream& f,  Vertex & v )
  { f >> (R2 &) v >> (Label &) v ; return f; }
  Vertex() : R2(),Label(){};
  Vertex(R2 P,int r=0): R2(P),Label(r){}
private:
  Vertex(const Vertex &);
  void operator=(const Vertex &);
};


// La classe BoundaryEdge (arêtes frontières)
class BoundaryEdge: public Label {
public:
  Vertex *vertices[2];
  void set(Vertex * v0,int i0,int i1,int r)
  { vertices[0]=v0+i0; vertices[1]=v0+i1; lab=r; }
  bool in(const Vertex * pv) const {return pv == vertices[0] || pv == vertices[1];}
  BoundaryEdge(){}; // constructor par défaut vide 
  void Draw() const;
  Vertex & operator[](int i) const {ASSERTION(i>=0 && i <2);
  return *vertices[i];}
  R length() const { R2 AB(*vertices[0],*vertices[1]);return sqrt((AB,AB));}
private:
  BoundaryEdge(const BoundaryEdge &);   // interdit la construction par copie
  void operator=(const BoundaryEdge &); // interdit l'affectation par copie 
};


// La classe Triangle (modélisation des triangles)
class Triangle: public Label {
  Vertex *vertices[3]; // tableau de trois pointeurs de type Vertex
public:
  R area;
  Triangle(){}; // constructeur par défaut vide
  Vertex & operator[](int i) const {
    ASSERTION(i>=0 && i <3);
    return *vertices[i];} // évaluation des pointeurs 
  void set(Vertex * v0,int i0,int i1,int i2,int r) {
    vertices[0]=v0+i0; vertices[1]=v0+i1; vertices[2]=v0+i2; 
    R2 AB(*vertices[0],*vertices[1]);
    R2 AC(*vertices[0],*vertices[2]);
    area = (AB^AC)*0.5;
    lab=r;
    ASSERTION(area>=0); }
  
  R2 Edge(int i) const {
    ASSERTION(i>=0 && i <3);
    return R2(*vertices[(i+1)%3],*vertices[(i+2)%3]);}// l'arête opposée au sommet i
  R2 H(int i) const { ASSERTION(i>=0 && i <3);
  R2 E=Edge(i);return E.perp()/(2*area);} // la hauteur
  R lenEdge(int i) const {
    ASSERTION(i>=0 && i <3);
    R2 E=Edge(i);return sqrt((E,E));}
private:
  Triangle(const Triangle &);       // interdit la construction par copie
  void operator=(const Triangle &); // interdit l'affectation  par copie
  
};



// La classe Mesh (modélisation du maillage)
class Mesh { public:
  int nt,nv,neb;
  R area;
  Vertex *vertices;
  Triangle *triangles;
  BoundaryEdge  *bedges;  
  
  Triangle & operator[](int i) const {return triangles[CheckT(i)];}
  Vertex & operator()(int i) const {return vertices[CheckV(i)];}
  inline Mesh(const char * filename); // lecture du fichier ".msh"
  int operator()(const Triangle & t) const {return CheckT(&t - triangles);}
  int operator()(const Triangle * t) const {return CheckT(t - triangles);}
  int operator()(const Vertex & v) const {return CheckV(&v - vertices);}
  int operator()(const Vertex * v) const{return CheckT(v - vertices);}
  int operator()(int it,int j) const {return (*this)(triangles[it][j]);}
  //  vérification des depassements de tableau 
  int CheckV(int i) const { ASSERTION(i>=0 && i < nv); return i;} 
  int CheckT(int i) const { ASSERTION(i>=0 && i < nt); return i;}

private:
  Mesh(const Mesh &);              // interdit la construction par copie
  void operator=(const Mesh &);    // interdit l'affectation  par copie
};

// Le constructeur de la classe Mesh
inline Mesh::Mesh(const char * filename)
{ // lecture du maillage
  int i,i0,i1,i2,ir;
  ifstream f(filename);
  if(!f) {cerr << "Mesh::Mesh Erreur a l'ouverture - fichier " << filename<<endl;exit(1);}
  cout << "Lecture du fichier  \"" <<filename<<"\""<<  endl;
  f >> nv >> nt >> neb ;
  cout << " Nb de sommets " << nv << " " << " Nb de triangles " 
       << nt << " Nb d'aretes frontiere " << neb << endl;
  assert(f.good() && nt && nv) ;
  triangles = new Triangle [nt];
  vertices = new Vertex[nv];
  bedges = new BoundaryEdge[neb];

  area=0;
  assert(triangles && vertices);
  for (i=0;i<nv;i++)    
    f >> vertices[i],assert(f.good());
  for (i=0;i<nt;i++) { 
    f >> i0 >> i1 >> i2 >> ir;
    assert(f.good() && i0>0 && i0<=nv && i1>0 && i1<=nv && i2>0 && i2<=nv);
    triangles[i].set(vertices,i0-1,i1-1,i2-1,ir); 
    area += triangles[i].area;}
  
  for (i=0;i<neb;i++) { 
    f >> i0 >> i1  >> ir;
    assert(f.good() && i0>0 && i0<=nv && i1>0 && i1<=nv );
    bedges[i].set(vertices,i0-1,i1-1,ir);}

   cout << " Fin lecture : aire du maillage = " << area <<endl;  
} 
