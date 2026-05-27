#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string>
#include <cmath>
#include <mpi.h>


class Comm{   // Communication internes
private:
  int* _NyLoc;
  int* _NxLoc;
  int _rank;
  MPI_Comm _comm_2D;
  int* _neighbor;
  int _nb_req;
  
public:
  Comm( int* NxLoc, int* NyLoc, MPI_Comm Comm_2D, int* neighbor): _NxLoc(NxLoc), _NyLoc(NyLoc), _comm_2D(Comm_2D),_neighbor(neighbor){
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);                        // Indicage des proc
  }

    
  MPI_Datatype DefMess(){     // On définit un nouveau type de données représentant une colonne de notre matrice 
    MPI_Datatype MPI_COL;
    MPI_Type_vector(_NyLoc[_rank],1,_NxLoc[_rank]+2,MPI_DOUBLE,&MPI_COL);
    MPI_Type_commit(&MPI_COL);
    return MPI_COL;
  }


  
  int AsynchroneComm(double* x, MPI_Request* request,MPI_Datatype MPI_COL) {  // Méthode de communication asynchrone entre les partitions juxatposées 
    int etiquette=100;
    int _nb_req = 0;

    for (int i=0;i < 8 ;i++) request[i]=MPI_REQUEST_NULL;


    MPI_Irecv(&x[0+1*(_NxLoc[_rank]+2)],1,MPI_COL,_neighbor[0],etiquette,_comm_2D,&request[_nb_req++]); 
    MPI_Irecv(&x[(_NxLoc[_rank]+1)+1*(_NxLoc[_rank]+2)],1,MPI_COL,_neighbor[1],etiquette,_comm_2D,&request[_nb_req++]);

    MPI_Isend(&x[1+1*(_NxLoc[_rank]+2)],1,MPI_COL,_neighbor[0],etiquette,_comm_2D,&request[_nb_req++]); // Envoi des valeurs 
    MPI_Isend(&x[_NxLoc[_rank]+1*(_NxLoc[_rank]+2)],1,MPI_COL,_neighbor[1],etiquette,_comm_2D,&request[_nb_req++]);   
    
    MPI_Irecv(&x[1+(_NyLoc[_rank]+1)*(_NxLoc[_rank] + 2) ],_NxLoc[_rank],MPI_DOUBLE,_neighbor[3],etiquette,_comm_2D,&request[_nb_req++]); 
    MPI_Irecv(&x[1+0*(_NxLoc[_rank] + 2)],_NxLoc[_rank],MPI_DOUBLE,_neighbor[2],etiquette,_comm_2D,&request[_nb_req++]);

    MPI_Isend(&x[1+_NyLoc[_rank]*(_NxLoc[_rank]+2)],_NxLoc[_rank],MPI_DOUBLE,_neighbor[3],etiquette,_comm_2D,&request[_nb_req++]); // Envoi des valeurs 
    MPI_Isend(&x[1+1*(_NxLoc[_rank]+2)],_NxLoc[_rank],MPI_DOUBLE,_neighbor[2],etiquette,_comm_2D,&request[_nb_req++]);

    return _nb_req;
  }

  void Wait(MPI_Request* request,int _nb_req){ // Attente de la fin des communications asynchrones    

    MPI_Waitall(_nb_req,request,MPI_STATUS_IGNORE); 
    
  }


  void SynchroneComm(double* x,MPI_Datatype MPI_COL){ // Méthode de communication synchrone entre les partitions juxatposées
    int etiquette=100;

    MPI_Sendrecv(&x[1+1*(_NxLoc[_rank]+2)],1,MPI_COL,_neighbor[0],etiquette,&x[(_NxLoc[_rank]+1)+1*(_NxLoc[_rank]+2)],1,MPI_COL,_neighbor[1],etiquette,_comm_2D,MPI_STATUS_IGNORE);

    MPI_Sendrecv(&x[_NxLoc[_rank]+1*(_NxLoc[_rank]+2)],1,MPI_COL,_neighbor[1],etiquette,&x[0+1*(_NxLoc[_rank]+2)],1,MPI_COL,_neighbor[0],etiquette,_comm_2D,MPI_STATUS_IGNORE);
     

    MPI_Sendrecv(&x[1+_NyLoc[_rank]*(_NxLoc[_rank]+2)],_NxLoc[_rank],MPI_DOUBLE,_neighbor[3],etiquette,&x[1],_NxLoc[_rank],MPI_DOUBLE,_neighbor[2],etiquette,_comm_2D,MPI_STATUS_IGNORE);
    
    // Échange des valeurs entre bornes juxtaposées

    MPI_Sendrecv(&x[1+1*(_NxLoc[_rank]+2)],_NxLoc[_rank],MPI_DOUBLE,_neighbor[2],etiquette,&x[1+(_NyLoc[_rank]+1)*(_NxLoc[_rank]+2)],_NxLoc[_rank],MPI_DOUBLE,_neighbor[3],etiquette,_comm_2D,MPI_STATUS_IGNORE);


  }

};


       
class Mesh{   // Classe réalisant notre méthode de calcul parallélisé 
private:
  double _Lx;                // Longueur du maillage en abscisse
  double _Ly;                // Longueur du maillage en ordonnée
  int _NbxElemtsGlob;        // Nombre de cellules en abscisse définies globalement
  int _nx;                   // Nombre réel d'élément sur une ligne 
  int _NbyElemtsGlob;        // Nombre de cellules en ordonnée définies globalement
  int _ny;                   // Nombre réel d'élement sur une colonne
  int* _debutx;              // Premier indice d'un proc sur l'axe des abscisses
  int* _debuty;              // Premier indice d'un proc sur l'axe des ordonnées
  int* _NbElemtsLoc;         // Vecteur contenant le nombre d'élements de chaque proc
  int* _NyLoc;               // Vecteur contenant le nombre de lignes de chaque proc
  int* _NxLoc;               // Vecteur contenant le nombre de colonnes de chaque proc
  int* _deplacement;         // Vecteur contenant le premier indice de chaque partition du maillage  
  int _nbp;                  // Nombre total de proc utilisés
  int _rank;                 // Indice permettant d'identifier chaque proc 
  double* _Init;             // Vecteur initial
  double* _Mesh0x;           // Maillage en abscisse initial
  double* _Mesh0y;           // Maillage en ordonnée initial 
  int* _dims;                // Nombre de proc par dimension
  int _neighbor[4];          // Rang des proc voisins
  MPI_Comm comm_2D;          // Communicateur dans un repère cartésien
  int* _coords;              // Coordonnées cartésiennes des prorcs
  int* _NbLocGhost;          // Nombre de cellules locales incluant les cellules fantômes
  
public:
  Mesh(const int& Nx,const int& Ny,const double& Lx,const double& Ly): _NbxElemtsGlob(Nx),_NbyElemtsGlob(Ny),_Lx(Lx), _Ly(Ly) {   // Constructeur de la classe
    MPI_Comm_size(MPI_COMM_WORLD,&_nbp);                         // Obtention du nombre de proc
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);                        // Indicage des proc
    _ny=_NbyElemtsGlob+1;
    _nx=_NbxElemtsGlob+1;
  }



  MPI_Comm DefComm(){  //On créer notre repère cartésien pour faciliter la communication
    int _ndims = 2;
    int _reorganisation;
    int _periode[_ndims];
    _coords = new int[2];
    _dims = new int[2];

    _dims[0]=0;
    _dims[1]=0;
    _periode[0]=false;
    _periode[1]=false;
    _reorganisation=false;

    MPI_Dims_create(_nbp, 2 ,_dims);
    if(_rank == 0) std::cout<<"dims[x]:"<<_dims[0]<<"\t dims[y]:"<<_dims[1]<<std::endl;; 
    
    MPI_Cart_create(MPI_COMM_WORLD,_ndims,_dims,_periode,_reorganisation,&comm_2D);
    MPI_Cart_coords(comm_2D,_rank,_ndims,_coords);

    return comm_2D;
  }


  void Defneighbor(){  // On identifie les voisins de notre proc
    int Left=0, Right=1, Bottom=2, Top=3;

    MPI_Cart_shift(comm_2D,0,1,&_neighbor[Left],&_neighbor[Right]);

    MPI_Cart_shift(comm_2D,1,1,&_neighbor[Bottom],&_neighbor[Top]);


  }

  int* Neighbor(){return _neighbor;}

  void computeNbElemtsLoc(){   // Définition du nombre de cellules sur chaque partition
    _NbElemtsLoc= new int[_nbp];
    _NyLoc= new int[_nbp];
    _NxLoc= new int[_nbp];

    
   
    _NyLoc[_rank] = _ny / _dims[1];
    if (_coords[1] < _ny % _dims[1] && _ny  % _dims[1] > 0) _NyLoc[_rank]++;
    _NxLoc[_rank] = _nx / _dims[0];
    if (_coords[0] < _nx % _dims[0] && _nx % _dims[0] > 0) _NxLoc[_rank]++;
    _NbElemtsLoc[_rank]=_NyLoc[_rank]*_NxLoc[_rank];

    MPI_Allgather(&_NbElemtsLoc[_rank],1,MPI_INT,_NbElemtsLoc,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&_NyLoc[_rank],1,MPI_INT,_NyLoc,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&_NxLoc[_rank],1,MPI_INT,_NxLoc,1,MPI_INT,MPI_COMM_WORLD);

  }
  int* NyLoc(){return _NyLoc;}
  int* NxLoc(){return _NxLoc;}

  void NbLocGhost(){             //Taille globale cellules fantômes comprises                 
    _NbLocGhost= new int[_nbp];
    _NbLocGhost[_rank] = (_NyLoc[_rank]+2)*(_NxLoc[_rank]+2) ;
  }
  
  void computeDeplacement(){  // Définiton des premiers indices de chaque cellules 
    _deplacement= new int[_nbp];
    _deplacement[0]=0;
    for(int j=1;j<_nbp;j++){
      _deplacement[j]=_deplacement[j-1]+_NbElemtsLoc[j-1];
    }
  }


  void setDebut(){
    _debuty = new int[_nbp];
    _debutx = new int[_nbp];
    
    _debuty[_rank] = _coords[1] * (_ny / _dims[1]) + std::min(_coords[1], _ny % _dims[1]);

    _debutx[_rank] = _coords[0] * (_nx / _dims[0]) + std::min(_coords[0], _nx % _dims[0]);
    MPI_Allgather(&_debuty[_rank],1,MPI_INT,_debuty,1,MPI_INT,MPI_COMM_WORLD);
	
    MPI_Allgather(&_debutx[_rank],1,MPI_INT,_debutx,1,MPI_INT,MPI_COMM_WORLD);
   
  }  // Définition du premier indice de chaque partition sans chaque maillage  

  int* Debuty() { return _debuty;}
  int* Debutx() { return _debutx;}


  double* AllocateMemory(){       
    return new double[_NbLocGhost[_rank]];   // Allocation mémoire d'une partition avec un point fantôme à chaque extrémité 
  }

  
  double* DefMeshx(){   // Défnition de chaque point du maillage global 
    _Mesh0x=new double[_nx];
    for ( int i=0; i<=_NbxElemtsGlob;i++) _Mesh0x[i]=i*_Lx/double(_NbxElemtsGlob);
    return _Mesh0x;
  }

  
  double* DefMeshy(){   // Défnition de chaque point du maillage global 
    _Mesh0y=new double[_ny];
    for ( int i=0; i<=_NbyElemtsGlob;i++) _Mesh0y[i]=i*_Ly/double(_NbyElemtsGlob);
    return _Mesh0y;
  }

  
  double* Initialize(Solution* Sol){  // Initialisatoion de nos partitions par la solution choisie 
    _Init= new double[_NbLocGhost[_rank]];
    for (int k =1; k<=_NyLoc[_rank];k++){
      for (int i=1; i<=_NxLoc[_rank];i++){
     	int I = i+k*(_NxLoc[_rank]+2);
	
     	_Init[I]=Sol->Initial(_Mesh0x[_debutx[_rank]+i-1],_Mesh0y[_debuty[_rank]+k-1]) ;
      }
    }
    
    
    return _Init;
  }

  double* AllocateCommMemory(){
    return new double[_NyLoc[_rank]];
  }

  

  void Next(double* Utemp, double* U, const double& dtloc){ // Calcul de notre partition au temps suivant
    const double dx = _Lx / double (_NbxElemtsGlob);
    const double dy = _Ly / double (_NbyElemtsGlob);
   
    for (int j = 2 ; j< _NyLoc[_rank];j++){
      for (int i = 2; i<_NxLoc[_rank]; i++){
	int I = i+j*(_NxLoc[_rank]+2);
	Utemp[I]= U[I] + dtloc / dx / dx * ( U[I+1]- 2 * U[I]+ U[I-1] );
	Utemp[I]+= dtloc / dy / dy * (U[ I+_NxLoc[_rank]+2 ]- 2* U[I] + U[ I-_NxLoc[_rank]-2 ]);
	// Application de la méthode d'Euler Explicite 
      }
    }
  }

  void BoundNext(double* Utemp, double* U, const double& dtloc){ // Calcul des nos bornes de chaque partition nécessitant les communications 
    const double dx = _Lx / double (_NbxElemtsGlob);
    const double dy = _Ly / double (_NbyElemtsGlob);

    
     for (int i = 1; i<= _NxLoc[_rank]; i++){
       int j = 1;   //Bottom 
       int I = i + j*(_NxLoc[_rank]+2);

       if (_neighbor[2]!=MPI_PROC_NULL){	 

	 if(i==1){
	   if (_neighbor[0]!=MPI_PROC_NULL ){
	   Utemp[I]=U[I] + dtloc /dx /dx * ( U[I+1] - 2*U[I] + U[I-1]);
	   Utemp[I]+= dtloc / dy / dy * (U[ I+ _NxLoc[_rank]+2 ]- 2* U[I] + U[ I- _NxLoc[_rank]-2 ]);
	   }
	   else Utemp[I]=U[I];
	 }
	 
	 else if(i==_NxLoc[_rank]){
	   if (_neighbor[1]!=MPI_PROC_NULL ){
	     Utemp[I]=U[I] + dtloc /dx /dx * ( U[I+1] - 2*U[I] + U[I-1]);
	     Utemp[I]+= dtloc / dy / dy * (U[ I+ _NxLoc[_rank]+2 ]- 2* U[I] + U[ I- _NxLoc[_rank]-2 ]);
	   }
	   else Utemp[I]=U[I];
	 }		
	 else {
	   Utemp[I]=U[I] + dtloc /dx /dx * ( U[I+1] - 2*U[I] + U[I-1]);
	   Utemp[I]+= dtloc / dy / dy * (U[ I+ _NxLoc[_rank]+2 ]- 2* U[I] + U[ I- _NxLoc[_rank]-2 ]);
	 }
       }
       else Utemp[I]= U[I];
     }


       
     for (int i = 1; i<= _NxLoc[_rank]; i++){
       int j = _NyLoc[_rank]; // Top
       int I= i+j*(_NxLoc[_rank]+2);
       
     if (_neighbor[3]!=MPI_PROC_NULL){	 

	 if(i==1){
	   if (_neighbor[0]!=MPI_PROC_NULL ){
	   Utemp[I]=U[I] + dtloc /dx /dx * ( U[I+1] - 2*U[I] + U[I-1]);
	   Utemp[I]+= dtloc / dy / dy * (U[ I+ _NxLoc[_rank]+2 ]- 2* U[I] + U[ I- _NxLoc[_rank]-2 ]);
	   }
	   else Utemp[I]=U[I];
	 }
	 
	 else if(i==_NxLoc[_rank]){
	   if (_neighbor[1]!=MPI_PROC_NULL ){
	     Utemp[I]=U[I] + dtloc /dx /dx * ( U[I+1] - 2*U[I] + U[I-1]);
	     Utemp[I]+= dtloc / dy / dy * (U[ I+ _NxLoc[_rank]+2 ]- 2* U[I] + U[ I- _NxLoc[_rank]-2 ]);
	   }
	   else Utemp[I]=U[I];
	 }		
	 else {
	   Utemp[I]=U[I] + dtloc /dx /dx * ( U[I+1] - 2*U[I] + U[I-1]);
	   Utemp[I]+= dtloc / dy / dy * (U[ I+ _NxLoc[_rank]+2 ]- 2* U[I] + U[ I- _NxLoc[_rank]-2 ]);
	 }
       }
       else Utemp[I]= U[I];
     }
     for (int j = 2; j< _NyLoc[_rank]; j++){
       int i = 1 ; //Left
       int I = i+j*(_NxLoc[_rank]+2);
      
      
       if (_neighbor[0]!=MPI_PROC_NULL){
	 Utemp[I]=U[I] + dtloc /dx /dx * ( U[I+1] - 2*U[I] + U[I-1]);
	 Utemp[I]+= dtloc / dy / dy * (U[ I+ _NxLoc[_rank]+2 ]- 2* U[I] + U[ I- _NxLoc[_rank]-2 ]);
       }
       
       else Utemp[I]=U[I];
     }
     for (int j = 2; j< _NyLoc[_rank]; j++){
       int i = _NxLoc[_rank]; // Right
       int I =  i+j*(_NxLoc[_rank]+2);
     	
       
       if (_neighbor[1]!=MPI_PROC_NULL){
	 Utemp[I]=U[I] + dtloc /dx /dx * ( U[I+1] - 2*U[I] + U[I-1]);
	 Utemp[I]+= dtloc / dy / dy * (U[ I+ _NxLoc[_rank]+2 ]- 2* U[I] + U[ I- _NxLoc[_rank]-2 ]);
	 


       }
       
       else Utemp[I]=U[I];
       
     }      
    
      
  }

  double* Copy(double* Utemp, double* U){  // Redéfiniton de U au temps t suivant
    for (int k =1; k<=_NyLoc[_rank];k++){
      for (int i=1; i<=_NxLoc[_rank];i++){
	int I = i+k*(_NxLoc[_rank]+2);
	U[I]=Utemp[I];
      }
    }
    return U;
  }



  void Write(const int& Ind,const double* const x, Solution* Sol,double& T){  // Écriture externe des solutions au temps T

    double* TempLoc=new double[_NxLoc[_rank]*_NyLoc[_rank]]; 
    for (int k =1; k<=_NyLoc[_rank];k++){
      for (int i=1; i<=_NxLoc[_rank];i++){
     	int I = i+k*(_NxLoc[_rank]+2);
	int Iloc = (i-1) + (k-1)*_NxLoc[_rank];
	TempLoc[Iloc]=x[I];
	//std::cout<<Iloc<<"\t"<<TempLoc[Iloc]<<std::endl;
      }
    }
    
    double* Xtmp = nullptr;
    double* X = nullptr;
    if(_rank == 0){
      X = new double[_nx*_ny];
      Xtmp = new double[_nx*_ny];
    }
    MPI_Gatherv(TempLoc,_NbElemtsLoc[_rank],MPI_DOUBLE,Xtmp,_NbElemtsLoc,_deplacement,MPI_DOUBLE,0,MPI_COMM_WORLD);  // Réunion de nos partition sur un unique proc
    delete[] TempLoc;
    
    if (_rank == 0){   
      
    for ( int r = 0; r < _nbp ; r++){
      for (int i= 0; i < _NxLoc[r]; i++){
	for (int j=0; j < _NyLoc[r];j++){
	  int I = _deplacement[r]+i + j*_NxLoc[r];
	  int J = _debutx[r] + i + _nx*(_debuty[r]+j);
	  X[J] = Xtmp[I];
	  //std::cout<<J<<"\t"<<Xtmp[I]<<std::endl;
	}
      }
    }
	  
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
  
