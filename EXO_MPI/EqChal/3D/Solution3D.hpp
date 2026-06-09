#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string>
#include <cmath>
#include <mpi.h>

class Solution{  
private:
  bool _Exist;  // Permet de dire si l'on connait la solution exacte ou non de la solution initiale associée

public:
  Solution(bool Exist) {  _Exist = Exist; }  // Constructeur de la classe

  virtual double Initial ( const double& x,const double& y,const double& z) {return 0 ; } // Fonction initiale héritant des classes dérivées


  virtual double Exact ( const double& x,const double& y,const double& z, const double& t){ 
    std::cerr<< "There is no exact function"<<std::endl;   // Fonction exacte héritant des classes dérivées
    return 0;
      } 
  bool HasExactSolution(){return _Exist;} // Vérification de l'existence de la solution exacte;
};
  

class Sin : public Solution {  // Classe dérivée de la classe Solution
private:
  int _k ;      // Diffusivité thermique
  double _Lx;   // Longueur du maillage en abscisse
  double _Ly;   // Longueur du maillage en ordonnée
  double _Lz;
public:
  Sin (const int& k,const double& Lx,const double& Ly,const double& Lz) : _k(k),_Lx(Lx),_Ly(Ly),_Lz(Lz), Solution(true){} // Constructeur de la classe

  double Initial(const double& x,const double& y,const double& z){return sin(_k*M_PI*x/_Lx)*sin(_k*M_PI*y/_Ly)*sin(_k*M_PI*z/_Lz);}   // Solution initiale

  double Exact(const double& x,const double& y,const double& z,const double& t) { return exp(-t* _k*_k* M_PI*M_PI * (1/_Lx/_Lx + 1/_Ly/_Ly + 1/_Lz/_Lz) )*Initial(x,y,z);} // Solution exacte
    
};


class Lin : public Solution { // Classe dérivée de la classe Solution

public:
  Lin ():Solution(true){}  // Constructeur de la classe

  double Initial(const double& x,const double& y,const double& z){return 12*(x+y+z)-17;} // Solution initiale

  double Exact(const double& x,const double& y,const double& z,const double& t) { return Initial(x,y,z);} //Solution exacte 	 
    
};

class Const : public Solution { // Classe dérivée de la classe Solution

public:
  Const ():Solution(true){}   // Constructeur de la classe

  double Initial(const double& x,const double& y,const double& z){return 4;} // Solution initiale

  double Exact(const double& x,const double& y,const double& z,const double& t) { return 4;} // Solution exacte 
    
};
