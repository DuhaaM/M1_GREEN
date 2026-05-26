#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string>
#include <cmath>
#include <mpi.h>
#include "Solution2D.hpp"
#include "Error2D.hpp"
#include "Mesh2D.hpp"


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

  Mesh.DefComm();
  Mesh.Defneighbor();
  Mesh.computeNbElemtsLoc(); // On calcule nos premières valeurs permettant nos calculs
  Mesh.NbLocGhost();
  Mesh.computeDeplacement();
  Mesh.setDebut();

  int* NxLoc=Mesh.NxLoc();
  int* NyLoc=Mesh.NyLoc();

  double* x = Mesh.DefMeshx();  // Définition du maillage global en abscisse
  double* y = Mesh.DefMeshy();  // Définition du maillage global en ordonnée

  double dx = Lx / double(Nx);// Pas d'espace
  double dy = Ly / double(Ny);
  double dtloc= 0.2 * std::min(dx,dy) * std::min(dx,dy);      // Pas de temps
  double t=0,time=0;
  int kT=0,Ind=0;
 
  double* Utemp= Mesh.AllocateMemory(); // Allocation mémoire du vecteur au temps t+1

 

  MPI_Request request[8];  // Création de nos identifiants de requête pour la méthode asynchrone
  
  
  double* U =  Mesh.Initialize(Sol);  // Initialisation du maillage en fonction de la solution
 
  MPI_Datatype MPI_COL = Mesh.DefMess();
       
       
  while( t < Tfinal ) {  // Boucle de calcul sur le temps 
  
    if ( t + dtloc > Tfinal) dtloc = Tfinal - t;   // Assurance de terminer au temps Tfinal 

    if (rang==0)temps_debut=MPI_Wtime();  // Marqueur de temps=0 pour le temps de calcul

    if ( method == 1) { // Utilisation de la méthode de communication asynchrone

      int nb_req = Mesh.AsynchroneComm(U,request,MPI_COL);                     // Réalistion des communications 
      Mesh.Next(Utemp,U,dtloc);                                                // Calcul au temps suivant
      Mesh.Wait(request,nb_req);                                               // Attente de la fin des communications 
      Mesh.BoundNext(Utemp,U,dtloc);                                           // Calcul des bords nécessitant les communications 
      Mesh.Copy(Utemp,U);                                                      // Passsage de U au temps suivant 

    }

    
    if ( method == 2) { // Utilisation de la méthode de communication synchrone

      Mesh.SynchroneComm(U,MPI_COL);              // Réalisation des communications 
      Mesh.Next(Utemp,U,dtloc);                                   // Calcul au temps suivant
      Mesh.BoundNext(Utemp,U,dtloc);                              // Calcul des bords nécessiant les communications
      Mesh.Copy(Utemp,U);                                         // Passage de U au temps suivant
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

  Error Error(Lx,Nx,Ny,K,Tfinal,x,y,Mesh.Debutx(),Mesh.Debuty(),NyLoc,NxLoc);

  
  
  if (rang == 0 ) std::cout<<"Computing time  = "<<time<<"s"<<std::endl;   //Sortie du temps de calcul  


  Error.ErrorInf(U,Sol);                   // Calcul de l'erreur en norme infinie 
  Error.ErrorL2(U,Sol);                    // Calcul de l'erreur en norme L2 

  

  
  MPI_Finalize();        // Fin de MPI
  return 0;

 
}
