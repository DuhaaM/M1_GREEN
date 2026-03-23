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

  virtual double Initial ( const double& x) {return 0 ; } // Fonction initiale héritant des classes dérivées


  virtual double Exact ( const double& x, const double& t){
    std::cerr<< "There is no exact function"<<std::endl;
    return 0;
      }  // Fonction exacte héritant des classes dérivées si ces dernières possèdent une solution exacte 

  bool HasExactSolution(){return _Exist;} // Vérification de l'existence de la solution exacte;
};
  

class Sin : public Solution {  // Classe dérivée de la classe Solution
private:
  int _k ;   // Diffusivité thermique
  double _L; // Longueur du maillage 

public:
  Sin (const int& k,const double& L) : _k(k),_L(L), Solution(true){} // Constructeur de la classe

  double Initial(const double& x){return sin(_k*M_PI*x/_L);}   // Solution initiale

  double Exact(const double& x,const double& t) { return exp(-t*_k*_k*M_PI*M_PI/_L/_L)*Initial(x);} // Solution exacte	 
    
};


class Lin : public Solution { // Classe dérivée de la classe Solution

public:
  Lin ():Solution(true){}  // Constructeur de la classe

  double Initial(const double& x){return 12*x-17;} // Solution initiale

  double Exact(const double& x,const double& t) { return Initial(x);} //Solution exacte 	 
    
};

class Const : public Solution { // Classe dérivée de la classe Solution

public:
  Const ():Solution(true){}   // Constructeur de la classe

  double Initial(const double& x){return 4;} // Solution initiale

  double Exact(const double& x,const double& t) { return 4;} // Solution exacte 
    
};

class Mesh{   // Classe réalisant notre méthode de calcul parallélisé 
private:
  double _L;                 //Longueur du maillage
  int _NbElemtsGlob;         // Nombre de cellules définies globalement 
  int _debut;                // Indice relatif au premier élément d'une partition par rapport au maillage initial
  int _fin;                  // Indice relatif au dernier élément d'une partition par rapport au maillage initial
  int* _NbElemtsLoc;         // Vecteur contenant le nombre d'élements de chaque proc
  int* _deplacement;         // Vecteur contenant le premier indice de chaque partition du maillage  
  int _nbp;                  // Nombre total de proc utilisés
  int _rank;                 // Indice permettant d'identifier chaque proc 
  double* _Init;             // Vecteur initial
  double* _Mesh0;            // Maillage initial 

public:
  Mesh(const int& N,const double& L): _NbElemtsGlob(N),_L(L) {   // Constructeur de la classe
    MPI_Comm_size(MPI_COMM_WORLD,&_nbp);                         // Obtention du nombre de proc
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);                        // Indicage des proc
  }


  void computeNbElemtsLoc(){   // Définition du nombre de cellules sur chaque partition
    _NbElemtsLoc= new int[_nbp];
    for(int j=0;j<_nbp;j++){
      _NbElemtsLoc[j]=(_NbElemtsGlob+1) / _nbp;
      if (j < (_NbElemtsGlob+1) % _nbp and (_NbElemtsGlob+1) % _nbp > 0) _NbElemtsLoc[j]++;
    }
  }

  
  void computeDeplacement(){  // Définiton des premiers indices de chaque cellules 
    _deplacement= new int[_nbp];
    _deplacement[0]=0;
    for(int j=1;j<_nbp;j++){
      _deplacement[j]=_deplacement[j-1]+_NbElemtsLoc[j-1];
    }
  }


  void setDebut(){ _debut= _deplacement[_rank];}  // Réecriture équivoque du premier inidice d'une partition 

  const int Debut() {return _debut;}
 
  void setFin(){ _fin = _debut+ _NbElemtsLoc[_rank]-1; } // Défnition du dernier indice d'une partition

  const int Fin() {return _fin;}


  double* AllocateMemory(){       
    return new double[_NbElemtsLoc[_rank]+2];   // Allocation mémoire d'une partition avec un point fantôme à chaque extrémité 
  }

  
  double* DefMesh(const double& L ){   // Défnition de chaque point du maillage global 
    _Mesh0=new double[_NbElemtsGlob+1];
    for ( int i=0; i<=_NbElemtsGlob;i++) _Mesh0[i]=i*_L/double(_NbElemtsGlob);
    return _Mesh0;
  }

  
  double* Initialize(Solution* Sol){  // Initialisatoion de nos partitions par la solution choisie 
    int j =1 ;
    _Init= new double[_NbElemtsLoc[_rank]+2];
      for (int i=_debut; i<=_fin;i++){
	_Init[j]=Sol->Initial(_Mesh0[i]) ;
	j++;
      }
      
      return _Init;
  }


  void AsynchroneComm(double* x, MPI_Request* request){  // Méthode de communication asynchrone entre les partitions juxatposées 
    int etiquette=100;

    for (int i=0;i < 4;i++) request[i]=MPI_REQUEST_NULL;

    
    if (_rank != _nbp-1) MPI_Irecv(&x[_NbElemtsLoc[_rank]+1],1,MPI_DOUBLE,_rank+1,etiquette,MPI_COMM_WORLD,&request[0]); // Définition d'une zone de réception
    if (_rank != 0)MPI_Irecv(&x[0],1,MPI_DOUBLE,_rank-1,etiquette,MPI_COMM_WORLD,&request[1]);

    if (_rank != _nbp-1) MPI_Isend(&x[_NbElemtsLoc[_rank]],1,MPI_DOUBLE,_rank+1,etiquette,MPI_COMM_WORLD,&request[2]); // Envoi des valeurs 
    if (_rank != 0) MPI_Isend(&x[1],1,MPI_DOUBLE,_rank-1,etiquette,MPI_COMM_WORLD,&request[3]);
  
  }

  void Wait(MPI_Request* request){  MPI_Waitall(4,request,MPI_STATUS_IGNORE);} // Attente de la fin des communications asynchrones  


  void SynchroneComm(double* x){ // Méthode de communication synchrone entre les partitions juxatposées
    int etiquette=100;

    if (_rank != _nbp-1) MPI_Sendrecv(&x[_NbElemtsLoc[_rank]],1,MPI_DOUBLE,_rank+1,etiquette,&x[_NbElemtsLoc[_rank]+1],1,MPI_DOUBLE,_rank+1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    // Échange des valeurs entre bornes juxtaposées

    if ( _rank != 0)  MPI_Sendrecv(&x[1],1,MPI_DOUBLE,_rank-1,etiquette,&x[0],1,MPI_DOUBLE,_rank-1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  }


  void Next(double* Utemp, double* U, const double& dtloc){ // Calcul de notre partition au temps suivant
    const double dx = _L / double (_NbElemtsGlob); 
    for (int i = 2; i<=_NbElemtsLoc[_rank]-1; i++){
      Utemp[i]= U[i] + dtloc / dx / dx * ( U[i+1]- 2 * U[i]+ U[i-1] ); // Application de la méthode d'Euler Explicite 
    }

    
  }

  void BoundNext(double* Utemp, double* U, const double& dtloc){ // Calcul des nos bornes de chaque partition nécessitant les communications 
    const double dx = _L / double (_NbElemtsGlob);

    
    if (_rank!=0) Utemp[1]=U[1] + dtloc /dx /dx * ( U[2] - 2*U[1] + U[0]);
    else Utemp[1]=U[1];

    if (_rank!=_nbp-1) Utemp[_NbElemtsLoc[_rank]] = U[_NbElemtsLoc[_rank]] + dtloc /dx /dx * ( U[_NbElemtsLoc[_rank]+1] -  2 * U[_NbElemtsLoc[_rank]] + U[_NbElemtsLoc[_rank]-1] );
    else Utemp[_NbElemtsLoc[_rank]]= U[_NbElemtsLoc[_rank]];
  }

  double* Copy(double* Utemp, double* U){  // Redéfiniton de U au temps t suivant
    for(int i = 1; i<=_NbElemtsLoc[_rank]; i++){
      U[i]=Utemp[i];
      }
    return U;
  }



  void Write(const int& Ind,const double* const x, Solution* Sol,double& T){  // Écriture externe des solutions au temps T 
    double* X = nullptr;
    if(_rank == 0)  X = new double[_NbElemtsGlob+1];
    MPI_Gatherv(&x[1],_NbElemtsLoc[_rank],MPI_DOUBLE,X,_NbElemtsLoc,_deplacement,MPI_DOUBLE,0,MPI_COMM_WORLD);  // Réunion de nos partition sur un unique proc
 
    if (_rank == 0){   
      std::stringstream mystream; 
      mystream<<"Data/Chal"<<Ind<<".csv"; // Définition du nom de fichier 
      
      std::ofstream File( mystream.str(), std::ofstream::out ); // Ouverture du fichier
      File << "#Indice , #Point maillage , #Solution obtenue , #Solution exacte  \n ";
      
      for(int i=0;i<=_NbElemtsGlob;i++){      
	File <<i<<","<< _Mesh0[i]<<","<<X[i]<<','<<Sol->Exact(_Mesh0[i],T)<<"\n";  //Écriture de l'indice / du point du maillage / de notre solution obtenue et de la solution exacte  
      }   
      File.close();  // Fermeture du fichier 
    }
  }


};
  


class Error{  //Classe calculant les erreurs de notre méthode
private:
  double _L;             // Longueur du maillage 
  int _NbElemtsGlob;     // Nombre d'éléments du maillage 
  int _k;                // Diffusivité thermique 
  double _Tfinal;        // Temps final 
  double* _Mesh0;        // Maillage défini
  int _fin;              // Dernier indice de chaque partition
  int _debut;            // Premier indice de chaque partition 
  int _rank;             // Indice de chaque proc 
  

public:
  Error(const double& L,const int& N, const int& k,const double& Tfinal,double* Init,const int& debut,const int& fin):_L(L), _NbElemtsGlob(N), _k(k), _Tfinal(Tfinal), _Mesh0(Init),_debut(debut),_fin(fin){
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);   //Constructeur de la classe
  } 


  void Verif(const int& kT,Solution* Sol){  // Propriété de la solution sinus vérifiant la convergence du schéma  
    double dx= _L /double(_NbElemtsGlob);
    double Normesqrt_U0=0,vp,Err;

    for (int i=0; i<=_NbElemtsGlob;i++){
      Normesqrt_U0 += dx*Sol->Initial(_Mesh0[i])*Sol->Initial(_Mesh0[i]);
    }
    
    vp= 4* (cos (_k*M_PI*dx) -1) / 5 +1;

    Err=abs( pow(vp,kT) - exp(-_Tfinal* _k*_k * M_PI*M_PI /_L/_L ) );

    if (_rank==0)std::cout<<"Erreur_cv="<<sqrt(Normesqrt_U0)*Err<<std::endl;

  }

  void ErrorInf(const double* Vect, Solution* Sol){  // Calcul de l'erreur en norme infinie
    double Errloc=0,Errglob=0; 
    int j=1;
    for (int i=_debut ; i<_fin ; i++){
      Errloc=std::max(Errloc,std::abs(Vect[j]-Sol->Exact(_Mesh0[i],_Tfinal)));
      j++;
    }
    MPI_Reduce(&Errloc,&Errglob,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD); //  Défnition du max entre tous les max de chaque proc 
    if(_rank==0) std::cout<<"ErreurInf="<<Errglob<<std::endl;
  }


  void ErrorL2(const double* Vect, Solution* Sol){ // Calcul approché de l'erreur L2 par une quadrature de Gauss en 2 points
    double Errloc=0,Errglob=0;
    double dx = _L/double(_NbElemtsGlob);
    double x_1=1 /sqrt(3), x_2=-1/sqrt(3);
    int j=1;
    for (int i=_debut; i < _fin ; i++){
      
      double xG1= (_Mesh0[i+1]+_Mesh0[i]) / 2 + (_Mesh0[i+1]-_Mesh0[i]) / 2 * x_1; // Défnition des points de Gauss
      double xG2= (_Mesh0[i+1]+_Mesh0[i]) / 2 + (_Mesh0[i+1]-_Mesh0[i]) / 2 * x_2;

      Errloc+= dx / 2 * (-Sol->Exact(xG1,_Tfinal) + (Vect[j+1]+Vect[j])/2 + (Vect[j+1]-Vect[j])/2 *x_1) * (-Sol->Exact(xG1,_Tfinal) + (Vect[j+1]+Vect[j])/2 + (Vect[j+1]-Vect[j])/2 *x_1);

      Errloc+= dx / 2 * (-Sol->Exact(xG2,_Tfinal) + (Vect[j+1]+Vect[j])/2 + (Vect[j+1]-Vect[j])/2 *x_2) * (-Sol->Exact(xG2,_Tfinal) + (Vect[j+1]+Vect[j])/2 + (Vect[j+1]-Vect[j])/2 *x_2);
      j++;
    }

    MPI_Reduce(&Errloc,&Errglob,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); // Somme des erreurs calculées sur chaque proc 
    
    if(_rank==0) std::cout<<"ErreurL2="<<sqrt(Errglob)<<std::endl;

  }
};
 
   
    
    



int main (int argc, char *argv[]){
  int N,K,kprint,method,Function,rang;
  double Tfinal,L;
  double temps_debut,temps_fin;

  // Définiton manuelle des variables 
  if (argc != 8) {
    std::cerr << "usage is : " << argv[0] << " <N=mesh size> <t=final time> <K=integer> <kprint=print step> <L=Upperbound> <Comm_Method : Non-Blocking = 1 or Sendrecv = 2 > <Function : 1 =  Sinus , 2 = Linear , 3 = Constant > "
	      << std::endl;
    return 0;
  }
  else {
  
    N=std::atoi(argv[1]); 
    Tfinal=std::atof(argv[2]);
    K=std::atoi(argv[3]);
    kprint=std::atoi(argv[4]);
    L=std::atof(argv[5]);
    method = std::atoi(argv[6]);
    Function  = std::atoi(argv[7]);
 }

  MPI_Init(&argc,&argv);    // Initialisation de MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rang); // On récupère les indices des proc pour écrire les sorties avec un seul proc

  Solution* Sol = nullptr;  // Création d'un pointeur nul vers notre solution 
 
  if (Function == 1) Sol = new Sin(K,L);

  if (Function == 2) Sol= new Lin() ;          //Définiton da la solution pointée

  if (Function == 3) Sol = new Const() ;
  
  Mesh Mesh(N,L);  

  Mesh.computeNbElemtsLoc(); // On calcule nos premières valeurs permettant nos calculs
  Mesh.computeDeplacement();
  Mesh.setDebut();
  Mesh.setFin();
 

  double* x = Mesh.DefMesh(L);  // Définition du maillage global

  double dx= L / double(N);      // Pas d'espace
  double dtloc= 0.4 * dx * dx;   // Pas de temps
  double t=0,time=0;
  int kT=0,Ind=0;
 
  double* Utemp= Mesh.AllocateMemory(); // Allocation mémoire du vecteur au temps t+1
 

  MPI_Request request[4];  // Création de nos identifiants de requête pour la méthode asynchrone
  
  
  double* U =  Mesh.Initialize(Sol);  // Initialisation du maillage en fonction de la solution 
  
   
   
  while( t < Tfinal ) {  // Boucle de calcul sur le temps 
  
    if ( t + dtloc > Tfinal) dtloc = Tfinal - t;   // Assurance de terminer au temps Tfinal 

    if (rang==0)temps_debut=MPI_Wtime();  // Marqueur de temps=0 pour le temps de calcul

    
    if ( method == 1) { // Utilisation de la méthode de communication asynchrone

      Mesh.AsynchroneComm(U,request);     // Réalistion des communications 
      Mesh.Next(Utemp,U,dtloc);           // Calcul au temps suivant
      Mesh.Wait(request);                 // Attente de la fin des communications 
      Mesh.BoundNext(Utemp,U,dtloc);      // Calcul des bords nécessitant les communications 
      Mesh.Copy(Utemp,U);                 // Passsage de U au temps suivant 

    }

    
    if ( method == 2) { // Utilisation de la méthode de communication synchrone

      Mesh.SynchroneComm(U);              // Réalisation des communications 
      Mesh.Next(Utemp,U,dtloc);           // Calcul au temps suivant
      Mesh.BoundNext(Utemp,U,dtloc);      // Calcul des bords nécessiant les communications
      Mesh.Copy(Utemp,U);                 // Passage de U au temps suivant
    }


    if(rang==0){
      temps_fin=MPI_Wtime();         // Marqueur de fin de temps pour le temps de calcul
      time+=temps_fin-temps_debut;   // Somme du temps de calcul à chaque passage dans la boucle 
    }


    

    

    if (kT % kprint == 0) {  // On écrit un fichier tous les kprint itérations 

      Mesh.Write(Ind,U,Sol,Tfinal); // Écriture du fichier 
      Ind++; 
     }
     
    t+=dtloc; // Passage au temps suivant 
    kT++;     // Nombre d'itérations 
  }

  Error Error(L,N,K,Tfinal,x,Mesh.Debut(),Mesh.Fin());

  
  
  if (rang == 0 ) std::cout<<"Computing time  = "<<time<<"s"<<std::endl;   //Sortie du temps de calcul  


  if (Function == 1)Error.Verif(kT,Sol);   // Vérification convergence de par les propriétés de la fonction sinus 
  Error.ErrorInf(U,Sol);                   // Calcul de l'erreur en norme infinie 
  Error.ErrorL2(U,Sol);                    // Calcul de l'erreur en norme L2 

  

  
  MPI_Finalize();        // Fin de MPI
  return 0;

 
}
  


  

  
