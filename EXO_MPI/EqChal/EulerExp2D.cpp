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

  virtual double Initial ( const double& x,const double& y) {return 0 ; } // Fonction initiale héritant des classes dérivées


  virtual double Exact ( const double& x,const double& y, const double& t){ 
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

public:
  Sin (const int& k,const double& Lx,const double& Ly) : _k(k),_Lx(Lx),_Ly(Ly), Solution(true){} // Constructeur de la classe

  double Initial(const double& x,const double& y){return sin(_k*M_PI*x/_Lx)*sin(_k*M_PI*y/_Ly);}   // Solution initiale

  double Exact(const double& x,const double& y,const double& t) { return exp(-t*_k*_k*M_PI*M_PI/_Lx/_Lx)*exp(-t*_k*_k*M_PI*M_PI/_Ly/_Ly)*Initial(x,y);} // Solution exacte	 
    
};


class Lin : public Solution { // Classe dérivée de la classe Solution

public:
  Lin ():Solution(true){}  // Constructeur de la classe

  double Initial(const double& x,const double& y){return 12*(x+y)-17;} // Solution initiale

  double Exact(const double& x,const double& y,const double& t) { return Initial(x,y);} //Solution exacte 	 
    
};

class Const : public Solution { // Classe dérivée de la classe Solution

public:
  Const ():Solution(true){}   // Constructeur de la classe

  double Initial(const double& x,const double& y){return 4;} // Solution initiale

  double Exact(const double& x,const double& y,const double& t) { return 4;} // Solution exacte 
    
};

class Mesh{   // Classe réalisant notre méthode de calcul parallélisé 
private:
  double _Lx;                // Longueur du maillage en abscisse
  double _Ly;                // Longueur du maillage en ordonnée
  int _NbxElemtsGlob;        // Nombre de cellules en abscisse définies globalement
  int _nx;                   // Nombre réel d'élémennt sur une ligne 
  int _NbyElemtsGlob;        // Nombre de cellules en ordonnée définies globalement
  int _debut;                // Indice relatif au premier élément d'une partition par rapport au maillage initial
  int _fin;                  // Indice relatif au dernier élément d'une partition par rapport au maillage initial
  int* _NbElemtsLoc;         // Vecteur contenant le nombre d'élements de chaque proc
  int* _NyLoc;               // Vecteur contenant le nombre de lignes de chaque proc
  int* _deplacement;         // Vecteur contenant le premier indice de chaque partition du maillage  
  int _nbp;                  // Nombre total de proc utilisés
  int _rank;                 // Indice permettant d'identifier chaque proc 
  double* _Init;             // Vecteur initial
  double* _Mesh0x;           // Maillage en abscisse initial
  double* _Mesh0y;           // Maillage en ordonnée initial 
  
public:
  Mesh(const int& Nx,const int& Ny,const double& Lx,const double& Ly): _NbxElemtsGlob(Nx),_NbyElemtsGlob(Ny),_Lx(Lx), _Ly(Ly) {   // Constructeur de la classe
    MPI_Comm_size(MPI_COMM_WORLD,&_nbp);                         // Obtention du nombre de proc
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);                        // Indicage des proc
    _nx=_NbxElemtsGlob+1;
  }


  void computeNbElemtsLoc(){   // Définition du nombre de cellules sur chaque partition
    _NbElemtsLoc= new int[_nbp];
    _NyLoc= new int[_nbp];
    
    for(int j=0;j<_nbp;j++){
      _NyLoc[j] = (_NbyElemtsGlob+1) / _nbp;
      if (j < (_NbyElemtsGlob+1) % _nbp and (_NbyElemtsGlob+1) % _nbp > 0) _NyLoc[j]++;
      _NbElemtsLoc[j]=_NyLoc[j]*(_nx);
    }
  }

  int* Nyloc(){return _NyLoc;}
  
  void computeDeplacement(){  // Définiton des premiers indices de chaque cellules 
    _deplacement= new int[_nbp];
    _deplacement[0]=0;
    for(int j=1;j<_nbp;j++){
      _deplacement[j]=_deplacement[j-1]+_NbElemtsLoc[j-1];
    }
  }


  void setDebut(){ _debut= _deplacement[_rank];}  // Réecriture équivoque du premier inidice d'une partition 

  const int Debut() {return _debut;}
 
  void setFin(){ _fin = _debut+ _NbElemtsLoc[_rank] -1; } // Défnition du dernier indice d'une partition

  const int Fin() {return _fin;}


  double* AllocateMemory(){       
    return new double[_NbElemtsLoc[_rank]+2*_nx];   // Allocation mémoire d'une partition avec un point fantôme à chaque extrémité 
  }

  
  double* DefMeshx(){   // Défnition de chaque point du maillage global 
    _Mesh0x=new double[_nx];
    for ( int i=0; i<=_NbxElemtsGlob;i++) _Mesh0x[i]=i*_Lx/double(_NbxElemtsGlob);
    return _Mesh0x;
  }

  
  double* DefMeshy(){   // Défnition de chaque point du maillage global 
    _Mesh0y=new double[_NbyElemtsGlob+1];
    for ( int i=0; i<=_NbyElemtsGlob;i++) _Mesh0y[i]=i*_Ly/double(_NbyElemtsGlob);
    return _Mesh0y;
  }
  
  double* Initialize(Solution* Sol){  // Initialisatoion de nos partitions par la solution choisie 
    int j =_nx ;
    _Init= new double[_NbElemtsLoc[_rank]+2*_nx];
    for (int k =0; k<_NyLoc[_rank];k++){
      for (int i=0; i<=_NbxElemtsGlob;i++){
	_Init[j++]=Sol->Initial(_Mesh0x[i],_Mesh0y[_debut/_nx+k]) ;
      }
    }
    
    return _Init;
  }


  void AsynchroneComm(double* x, MPI_Request* request){  // Méthode de communication asynchrone entre les partitions juxatposées 
    int etiquette=100;

    for (int i=0;i < 4;i++) request[i]=MPI_REQUEST_NULL;

    
    if (_rank != _nbp-1) MPI_Irecv(&x[_NbElemtsLoc[_rank]+_nx],_nx,MPI_DOUBLE,_rank+1,etiquette,MPI_COMM_WORLD,&request[0]); // Définition d'une zone de réception
    if (_rank != 0)MPI_Irecv(&x[0],_nx,MPI_DOUBLE,_rank-1,etiquette,MPI_COMM_WORLD,&request[1]);

    if (_rank != _nbp-1) MPI_Isend(&x[_NbElemtsLoc[_rank]],_nx,MPI_DOUBLE,_rank+1,etiquette,MPI_COMM_WORLD,&request[2]); // Envoi des valeurs 
    if (_rank != 0) MPI_Isend(&x[_nx],_nx,MPI_DOUBLE,_rank-1,etiquette,MPI_COMM_WORLD,&request[3]);
  
  }

  void Wait(MPI_Request* request){  MPI_Waitall(4,request,MPI_STATUS_IGNORE);} // Attente de la fin des communications asynchrones  


  void SynchroneComm(double* x){ // Méthode de communication synchrone entre les partitions juxatposées
    int etiquette=100;

    if (_rank != _nbp-1) MPI_Sendrecv(&x[_NbElemtsLoc[_rank]],_nx,MPI_DOUBLE,_rank+1,etiquette,&x[_NbElemtsLoc[_rank]+_nx],_nx,MPI_DOUBLE,_rank+1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    // Échange des valeurs entre bornes juxtaposées

    if ( _rank != 0)  MPI_Sendrecv(&x[_nx],_nx,MPI_DOUBLE,_rank-1,etiquette,&x[0],_nx,MPI_DOUBLE,_rank-1,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  }


  void Next(double* Utemp, double* U, const double& dtloc){ // Calcul de notre partition au temps suivant
    const double dx = _Lx / double (_NbxElemtsGlob);
    const double dy = _Ly / double (_NbyElemtsGlob);
    for (int j =1 ; j< _NyLoc[_rank];j++){
      for (int i = 1; i<_NbxElemtsGlob; i++){
	int I = i+j*_nx;
	Utemp[I]= U[I] + dtloc / dx / dx * ( U[I+1]- 2 * U[I]+ U[I-1] );
	Utemp[I]+= dtloc / dy / dy * (U[ I+_nx ]- 2* U[I] + U[ I-_nx ]);
	// Application de la méthode d'Euler Explicite 
      }
    } 
  }

  void BoundNext(double* Utemp, double* U, const double& dtloc){ // Calcul des nos bornes de chaque partition nécessitant les communications 
    const double dx = _Lx / double (_NbxElemtsGlob);
    const double dy = _Ly / double (_NbyElemtsGlob);
    
    for (int i = 1; i<_NbxElemtsGlob; i++){
      int I = i + _nx;

      if (_rank!=0){
	Utemp[I]=U[I] + dtloc /dx /dx * ( U[I+1] - 2*U[I] + U[I-1]);
	Utemp[I]+= dtloc / dy / dy * (U[ I+_nx ]- 2* U[I] + U[ I-_nx ]);
      }
      else Utemp[I]=U[I];

      I= i + _NyLoc[_rank]*_nx;

      if (_rank!=_nbp-1){
	Utemp[I] = U[I] + dtloc /dx /dx * ( U[I+1] -  2 * U[I] + U[I-1] );
	Utemp[I]+= dtloc / dy / dy * (U[ I+_nx ]- 2* U[I] + U[ I-_nx ]);
      }
      
      else Utemp[I]= U[I];
    }
  }

  double* Copy(double* Utemp, double* U){  // Redéfiniton de U au temps t suivant
    for (int j = 1 ; j<= _NyLoc[_rank];j++){
      for (int i = 1; i<_NbxElemtsGlob; i++){
	int I = i+j* _nx;
	U[I]=Utemp[I];
      }
    }
    return U;
  }



  void Write(const int& Ind,const double* const x, Solution* Sol,double& T){  // Écriture externe des solutions au temps T 
    double* X = nullptr;
    if(_rank == 0)  X = new double[_nx*(_NbyElemtsGlob+1)];
    MPI_Gatherv(&x[_nx],_NbElemtsLoc[_rank],MPI_DOUBLE,X,_NbElemtsLoc,_deplacement,MPI_DOUBLE,0,MPI_COMM_WORLD);  // Réunion de nos partition sur un unique proc
 
    if (_rank == 0){   
      std::stringstream mystream; 
      mystream<<"Data2D/Chal"<<Ind<<".csv"; // Définition du nom de fichier 
      
      std::ofstream File( mystream.str(), std::ofstream::out ); // Ouverture du fichier
      File << "#Indice,#Point_maillage1,#Point_maillage2,#Solution_obtenue,#Solution_exacte"<<"\n";

      int k = 0;

      for(int j=0;j<=_NbyElemtsGlob;j++){
	for( int i=0; i<=_NbxElemtsGlob;i++){
	  File <<k<<","<< _Mesh0x[i]<<","<<_Mesh0y[j]<<","<<X[k]<<','<<Sol->Exact(_Mesh0x[i],_Mesh0y[j],T)<<"\n";  //Écriture de l'indice / du point du maillage / de notre solution obtenue et de la solution exacte
	  k++;
	}
      }

      File.close();  // Fermeture du fichier 

    }
  }


};
  


class Error{  //Classe calculant les erreurs de notre méthode
private:
  double _L;             // Longueur du maillage 
  int _NbxElemtsGlob;    // Nombre d'éléments du maillageœ
  int _NbyElemtsGlob;    // Nombre d'éléments du maillage 
  int _k;                // Diffusivité thermique 
  double _Tfinal;        // Temps final 
  double* _Mesh0x;       // Maillage définiœ
  double* _Mesh0y;       // Maillage défini
  int _fin;              // Dernier indice de chaque partition
  int _debut;            // Premier indice de chaque partition 
  int _rank;             // Indice de chaque proc
  int* _NyLoc;

public:
  Error(const double& L,const int& Nx,const int& Ny, const int& k,const double& Tfinal,double* Initx ,double* Inity ,const int& debut,const int& fin,int* Nyloc):_L(L), _NbxElemtsGlob(Nx), _NbyElemtsGlob(Ny), _k(k), _Tfinal(Tfinal), _Mesh0x(Initx),_Mesh0y(Inity),_debut(debut),_fin(fin),_NyLoc(Nyloc){
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);   //Constructeur de la classe
  } 


  void Verif(const int& kT,Solution* Sol){  // Propriété de la solution sinus vérifiant la convergence du schéma  
    double dx= _L /double(_NbxElemtsGlob);
    double Normesqrt_U0=0,vp,Err;

    for (int i=0; i<=_NbxElemtsGlob;i++){
      Normesqrt_U0 += dx*Sol->Initial(_Mesh0x[i],0)*Sol->Initial(_Mesh0x[i],0);
    }
    
    vp= 4* (cos (_k*M_PI*dx) -1) / 5 +1;

    Err=abs( pow(vp,kT) - exp(-_Tfinal* _k*_k * M_PI*M_PI /_L/_L ) );

    if (_rank==0)std::cout<<"Erreur_cv="<<sqrt(Normesqrt_U0)*Err<<std::endl;

  }

  void ErrorInf(const double* Vect, Solution* Sol){  // Calcul de l'erreur en norme infinie
    double Errloc=0,Errglob=0; 
    int j=_NbxElemtsGlob+1;
    for (int k =0; k<_NyLoc[_rank];k++){
      for (int i=0; i<=_NbxElemtsGlob;i++){
	Errloc=std::max(Errloc,std::abs(Vect[j]-Sol->Exact(_Mesh0x[i],_Mesh0y[_debut/(_NbxElemtsGlob+1)+k],_Tfinal)));
	j++;
      }
    }
    MPI_Reduce(&Errloc,&Errglob,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD); //  Défnition du max entre tous les max de chaque proc 
    if(_rank==0) std::cout<<"ErreurInf="<<Errglob<<std::endl;
  }


  void ErrorL2(const double* Vect, Solution* Sol){ // Calcul approché de l'erreur L2 par une quadrature de Gauss en 2 points
    double Errloc=0,Errglob=0;
    double dx = _L/double(_NbxElemtsGlob);
    double x_1=1 /sqrt(3), x_2=-1/sqrt(3);
    int j=1;
    for (int i=_debut; i < _fin ; i++){
      
      double xG1= (_Mesh0x[i+1]+_Mesh0x[i]) / 2 + (_Mesh0x[i+1]-_Mesh0x[i]) / 2 * x_1; // Défnition des points de Gauss
      double xG2= (_Mesh0x[i+1]+_Mesh0x[i]) / 2 + (_Mesh0x[i+1]-_Mesh0x[i]) / 2 * x_2;

      Errloc+= dx / 2 * (-Sol->Exact(xG1,0,_Tfinal) + (Vect[j+1]+Vect[j])/2 + (Vect[j+1]-Vect[j])/2 *x_1) * (-Sol->Exact(xG1,0,_Tfinal) + (Vect[j+1]+Vect[j])/2 + (Vect[j+1]-Vect[j])/2 *x_1);

      Errloc+= dx / 2 * (-Sol->Exact(xG2,0,_Tfinal) + (Vect[j+1]+Vect[j])/2 + (Vect[j+1]-Vect[j])/2 *x_2) * (-Sol->Exact(xG2,0,_Tfinal) + (Vect[j+1]+Vect[j])/2 + (Vect[j+1]-Vect[j])/2 *x_2);
      j++;
    }

    MPI_Reduce(&Errloc,&Errglob,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); // Somme des erreurs calculées sur chaque proc 
    
    if(_rank==0) std::cout<<"ErreurL2="<<sqrt(Errglob)<<std::endl;

  }
};
 
   
    
    



int main (int argc, char *argv[]){
  int Nx,Ny,K,kprint,method,Function,rang;
  double Tfinal,Lx,Ly;
  double temps_debut,temps_fin;



  MPI_Init(&argc,&argv);    // Initialisation de MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rang); // On récupère les indices des proc pour écrire les sorties avec un seul proc

  // Définiton manuelle des variables 
  if (argc != 10) {
    std::cerr << "usage is : " << argv[0] << " <Nx=mesh size> <Ny=mesh size> <t=final time> <K=integer> <kprint=print step> <Lx=Upperbound> <Ly=Upperbound> <Comm_Method : Non-Blocking = 1 or Sendrecv = 2 > <Function : 1 =  Sinus , 2 = Linear , 3 = Constant > "
	      << std::endl;
    return 0;
  }
  else {
  
    Nx=std::atoi(argv[1]);
    Ny=std::atoi(argv[2]);
    Tfinal=std::atof(argv[3]);
    K=std::atoi(argv[4]);
    kprint=std::atoi(argv[5]);
    Lx=std::atof(argv[6]);
    Ly=std::atof(argv[7]);
    method = std::atoi(argv[8]);
    Function  = std::atoi(argv[9]);
 }

  Solution* Sol = nullptr;  // Création d'un pointeur nul vers notre solution 
 
  if (Function == 1) Sol = new Sin(K,Lx,Ly);

  if (Function == 2) Sol= new Lin() ;          //Définiton da la solution pointée

  if (Function == 3) Sol = new Const() ;
  
  Mesh Mesh(Nx,Ny,Lx,Ly);  

  Mesh.computeNbElemtsLoc(); // On calcule nos premières valeurs permettant nos calculs
  Mesh.computeDeplacement();
  Mesh.setDebut();
  Mesh.setFin();
 

  double* x = Mesh.DefMeshx();  // Définition du maillage global en abscisse
  double* y = Mesh.DefMeshy();  // Définition du maillage global en ordonnée

  double dx = Lx / double(Nx);// Pas d'espace
  double dy = Ly / double(Ny);
  double dtloc= 0.2 * std::min(dx,dy) * std::min(dx,dy);      // Pas de temps
  double t=0,time=0;
  int kT=0,Ind=0;
 
  double* Utemp= Mesh.AllocateMemory(); // Allocation mémoire du vecteur au temps t+1
 

  MPI_Request request[4];  // Création de nos identifiants de requête pour la méthode asynchrone
  
  
  double* U =  Mesh.Initialize(Sol);  // Initialisation du maillage en fonction de la solution

  if( rang==1 ){
    for (int i = 0; i<20;i++) std::cout<<i<<"\t"<<U[i]<<std::endl;
  }
  
   
   
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

      Mesh.Write(Ind,U,Sol,t); // Écriture du fichier 
      Ind++; 
     }
     
    t+=dtloc; // Passage au temps suivant 
    kT++;     // Nombre d'itérations 
  }

  Error Error(Lx,Nx,Ny,K,Tfinal,x,y,Mesh.Debut(),Mesh.Fin(),Mesh.Nyloc());

  
  
  if (rang == 0 ) std::cout<<"Computing time  = "<<time<<"s"<<std::endl;   //Sortie du temps de calcul  


  if (Function == 1)Error.Verif(kT,Sol);   // Vérification convergence de par les propriétés de la fonction sinus 
  Error.ErrorInf(U,Sol);                   // Calcul de l'erreur en norme infinie 
  Error.ErrorL2(U,Sol);                    // Calcul de l'erreur en norme L2 

  

  
  MPI_Finalize();        // Fin de MPI
  return 0;

 
}
  


  

  
