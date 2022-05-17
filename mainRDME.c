// 
// Monte Carlo 3D simulations of mesoscopic reaction-diffusion kinetics 
// The algorithm was implemented for 3D 
// Model for A=antibody V=virion to produce AV, A=Fab-free  V=LUVs+0Fab  
// A+V   <==> A1V
// A+A1V <==> A2V
// .
// .
// A+AmV <==> AmV
//
// Compile with: make
// Usage: ./rdme
// last modification : may 1 2022 LV lupevillegaslopez@gmail.com
// The program outputs 3 files with the positions of all the particles at each event time 
// "particlestimeEvolution.txt":           total number of particles of each specie:  t-A-V-A1V-A2V-...AmV
// "timeParticlest.txt":                   the particle in a diffusion event at time t  : "particle-t-x-y-z-especie"
// "timeEvolutionC.txt":                   the particle in a reaction event at time t  : "t-x-y-z-V-A-VA-reaction"
//
//Libreries  
#include <iostream>  	//cout, printf
#include <iomanip>
#include <vector>    	//vector class
#include <fstream>   	//ofstream, ifstream for file IO
#include <cmath>     	//sqrt for rdf
#include <algorithm> 	//sort algorithm for heal_list
#include <stdio.h>    	//printf 
#include <math.h>     	//exp, floor, log , sqrt
#include <stdlib.h>     //rand
#include <random>       //Gaussian noise
#include <time.h>       //clock_t, 
#include "ran2.h"      //Random numbers. Uniform distribution
using namespace std;   //funciones std.names


//Functions
int init();                                         //Initialize  
double Myrand();                                    //random numbers
void make_connectivity3D();                         //create the connectivity matrix 3D
void make_configuration(int i, int icell);          //create the configuration matrix
void make_linked_list3D();                          //Fill head and list with the current system state
void make_rateMatrix();                             //Create rate matrix
void make_rateMatrix3D();                           //Create rate matrix for 3 species
void make_event_q();                                //Create event queue array. First event time 
void make_q_matrix();                               //Make Q-array
void reaction_event(int sv, int sw);                //Reaction event
void diffusion_event(int sv, int sw);               //Diffusion event
void sort(int num_p);                               //sort array of 2xn variables
void apply_pbc3D(double *r);                        //Periodic Boundary Conditions initial positions 
double change_rate3D(double k, double D);           //Change rate constant (3D) k-->q Mesoscopic constants
int getcell3D(double *pos, int i);                 //Compute cell coordinates in cell, return cell index, save new positions 
void get_nextEvent(int sv);                        //Compute next event and update event matrix
void get_nextEventD(int nsv, int sv, double rs_old, double t_alpha);  //Compute next event and update event matrix
void update_rateMatrix3D(int svN);                         //Update rate matrix
void updateEvent();                                //Update event
void make_positionCell();                          //Make array positionCell
void move_particles(int nsv, int particle);        //save new positions of molecules
int get_molecule(int sv, int especie);             //get the molecule that moves
int get_Reaction(int sv, double q1, double q2);    //get the number of reaction
int get_Diffusion(int sv);                         //Get the particle that diffuses
void update_headList();                             //Update head and list
void update_moleculesR(int sv, int r, int mol1, int mol2, int mol3, int sav, int sw, int reaction); //Update mol in reactions
void write_particlestimeEvolution(const char* filename, double t);              //write total number of molecules
void write_particles(const char* filename, double t, int i);                    //Write particles in an event
void write_allparticles(const char* filename, double t);                        //Write all particles
void write_timeEvolutionC(const char* filename, double t);                      //write positions of AV
void writeEvent();                                  //Write event
void writeMolecules();                              //Write molecules
void writeSpecies();                                //Write species
void writeConfiguration();                          //Write configuration
void writeSummary();                                //Write summary 
void readEvent();                                   //Read event
void readMolecules();                               //Read molecules
void readSpecies();                                 //Read species
void readConfiguration();                           //Read configuration
void readParameters();                              //Read initial parameters
void writeInfo();                                   //Write info
void readInfo();                                    //Read info
void writeMolC();                                   //Write molecules C
void readMolC();                                    //Read molecules C


//***Parameters********
int    M;                                           //Number of species
int    dimension = 3;                               //2:2D dimension, 3: 3D dimension
int    N_A;                                         //total Number of particles of type Fab-free
int    N_V;                                         //total Number of particles of type LUV-free
int    N;                                           //total Number of particles
int    N_C = 9000000;                               //total number of reactions
double L;                                           // simulation volume size (3umx3umx3um)
int    NCeldas1D;                                   // Number of subvolumenes 
int    NCeldas2D;                                   // Number of subvolumenes 
int    Nsv;                                         // Total number of subvolumes
double l ;                                          // subvolume size (Ex: size(l)=0.5um, L=5um, NCeldas1D = 10)
double Av = 6.02214076e23;                          // Avogadro constant
double Molec2Molar = Av/1e15;                       // Change Molecules/um^3  to Molar units
double Molar2Molec = 1e15/Av;                       // Change Molar to Molecules/um^3 units
double pvol;                                        // Volumen-1   1/um^3 
int    R;                                           // number of reactions
double k1;                                          // s-1M-1   rate  V+A --> VA  ; kon [s-1M-1]
double k2;                                          // s-1     rate   VA --> V+ A  ;koff[s-1]
double D_A;                                         // Diffusion constant of A (um2/s)
double D_V;                                         // Diffusion constant of V, A1V, A2V, A3V...AmV (um2/s)
double rho_A;                                       //radius of Fab (um)
double rho_V;                                       //radius of LUV (um)
double rho;                                         //reaction radius
int nsteps;                                        //Number steps to perform
int save_freq;                                      //Write results every save_freq steps
int r1;                                             //number of reactions A+B-->C
int r2;                                             //number of reactions C-->A+B
int r12;                                            //current number of reaction
double tclock;                          
long seed;                                          //seed

vector<double> particles;        //Particle positions, stored aligned as x1,y1,x2,y2,x3,y3,..
vector<int> molecules;           //new positions of particles according to cells
vector<int> connectivity;        //Connectivity matrix, stored aligned as a11,a12,a13,a14,a15,a16,a21,a22...
vector<int> configuration;       //Configuration matrix, stored number of particles
vector<int> list, head, species; //List, head, and species vector
vector<double> rate;             //rate matrix
vector<double> event_queue;      //event queue array 
vector<double> event;            //event queue array (subvolumes estan ordenados en forma ascendente)
vector<double> parameters;       //List of parameters 
vector<int> q_array;             //Q-array
vector<int> positionCell;        //position of the cell
vector<int> moleculesC;          //positions of particles C 



//***********************MAIN*****************

int main(){
    
    clock_t treloj;
    treloj = clock();
    
    seed = (long)time(NULL);   
    // Initialize the random generator ran2 with NEGATIVE argument 
    ran2(-seed);        // This must be done only once! 
    
    int iter, sv, sw;
    double t, t0,  next, aleatorio;
    t0 = 0, r1 = 0, r2 = 0, r12 =0;
    
    // 1 - Initialize (Make configuration matrix)
    sw = init();
    if(sw==0){
        // 2. Make connectivity matrix
        make_connectivity3D();
        make_positionCell();
        // 3. Make rate_queue matrix
        make_rateMatrix3D();
        // 4. Make event queue
        make_event_q();
        // 5. Make q-array
        make_q_matrix();
        write_allparticles("timeParticlest.txt", t0);
        write_particlestimeEvolution("particlestimeEvolution.txt", t0);
    }

    if(sw==1){
        make_connectivity3D();
        make_positionCell();
        make_rateMatrix3D();
        t0 = event[1];
        write_allparticles("timeParticlest1.txt", t0);
        write_particlestimeEvolution("particlestimeEvolution1.txt", t0);  
    }
    
    // 6. Next event
    for(iter=1; iter<=nsteps; iter++)
    {             
        //time 
        t = event[1];
        //subvolume=lambda 
        sv = (int) (event[0]);  
        
        //Probability for next event
        next = rate[3*sv]/rate[3*sv+2];
        aleatorio = Myrand();        

        //Reaction event
        if(aleatorio < next)
        {
            //printf("reaction: %g %d \n",t, iter);
            reaction_event(sv, sw);
            if(sw==0)
            {
            write_particlestimeEvolution("particlestimeEvolution.txt", t);
			write_timeEvolutionC("timeEvolutionC.txt", t);
            }
            
            if(sw==1)
            {
            write_particlestimeEvolution("particlestimeEvolution1.txt", t);
			write_timeEvolutionC("timeEvolutionC1.txt", t);
            }
            //next event
            t = event[1];
        }
        //Diffusion event
        else
        {
            diffusion_event(sv, sw);
            //next event
            t = event[1]; 
        }
      
        if(iter%(save_freq)==0)
        {
            writeEvent();
            writeMolecules();	
            writeSpecies();
            writeConfiguration();
            writeInfo();
            writeMolC();
        }
        
    }
    treloj = clock() - treloj;
    tclock = ((float)treloj)/CLOCKS_PER_SEC;
    writeInfo();
    writeSummary();
    return 0;
}
  
//Initialize the variables and generate random positions
int init(){
    
    //Get initial parameters
    parameters.resize(13,0.0);
    readParameters();
    N_A = parameters[0];
    N_V = parameters[1];
    D_A = parameters[2];                        
    D_V = parameters[3];   
    rho_A = parameters[4]/1000;                 
    rho_V = parameters[5]/1000;  
    M = parameters[6] + 2;  
    L = parameters[7];
    NCeldas1D = parameters[8];            
    k1 = parameters[9];
    k2 = parameters[10];
    nsteps    = parameters[11];         
    save_freq = parameters[12];N_A = parameters[0];
    N_V = parameters[1];
    D_A = parameters[2];                        
    D_V = parameters[3];   
    rho_A = parameters[4]/1000;                 
    rho_V = parameters[5]/1000;  
    M = parameters[6] + 2;  
    L = parameters[7];
    NCeldas1D = parameters[8];            
    k1 = parameters[9];
    k2 = parameters[10];
    nsteps    = parameters[11];         
    save_freq = parameters[12];
    N = N_A + N_V;                
    rho = rho_A + rho_V; 
    NCeldas2D = NCeldas1D*NCeldas1D;     
    Nsv = NCeldas2D*NCeldas1D;                                                  
    l = L/NCeldas1D;                            
    pvol= 1/(l*l*l);                           
    R = 2*(M-2);                               


    //initialize
    particles.resize(dimension*N);               //vector of positions 
    molecules.resize(dimension*(N+N_C+1), 0);    //vector of new positions
    moleculesC.resize(6*(N_C+1),0);              //vector molecules C
    connectivity.resize(2*dimension*Nsv);        //connectivity matrix
    configuration.resize(Nsv*M,0);               //configuration matrix
    species.resize(N+N_C+1,M+1);                 //species vector
    rate.resize(3*Nsv,0.0);                     //rate matrix
    event_queue.resize(2*Nsv,0.0);              //event queue array   
    event.resize(2*Nsv,0.0);                    //event array
    q_array.resize(Nsv,0);                      //q-array
    positionCell.resize((dimension+1)*Nsv,0);   //array positiones of the subvolumes
    

    // We store in head (from 1 to ncells=Nsv) and list (from 1 to N).
    head.clear();
    list.clear();
    head.resize(Nsv+1, 0);
    list.resize(N+N_C+1, 0);
    FILE *f;
	f = fopen("event.txt", "r");
    int sw;
	int Nlast;
    
    if(f == NULL){
        
        //Initial position. Particles are localized in a box.
        for(int i=0; i<N; i++)
        {
            particles[3*i]   = Myrand()*L;
            particles[3*i+1] = Myrand()*L;
            particles[3*i+2] = Myrand()*L;
        }
    
		//Distribution 
        // A = 0 ;  A0V= 1; A1V=2; A2V=3; 
        for(int i=1; i<N_A+1; i++)
        {
            species[i]   = 0; 
        }
		Nlast = N_A;
        
		for(int i=Nlast+1; i<=Nlast+N_V;  i++)
        {
            species[i]   = 1; 
        }
		
        //Save in head and list
        make_linked_list3D();
        write_allparticles("initial.txt", 0);
        sw = 0;
    }
    else
    {
        readInfo();
        readEvent();
        updateEvent();
        readMolecules();
        readConfiguration();
        readMolC();
        readSpecies();
        update_headList();
        sw = 1;
    }
    
   return sw; 
}


//Make Connectivity Matrix 3D
void make_connectivity3D(){
    int dim_connec = 2*dimension*Nsv;
    int aux1 = NCeldas2D-NCeldas1D;
    int aux2 = NCeldas2D;
    int aux3 = 0;
    int aux4 = NCeldas1D;
    int aux5 = NCeldas2D;    
    int aux6 = NCeldas1D-1;
    int aux7 = NCeldas2D-NCeldas1D;

    vector<int> connec;  
    connec.resize(dim_connec,0);   

    for(int fila=0;fila<Nsv; fila++)
    {
      if( fila%NCeldas1D == NCeldas1D-1 )
      {
        connec[ 6*fila + 0 ] = fila -aux6;
      }
      else
      {
        connec[ 6*fila + 0 ] = fila + 1;
      }
      if( fila%NCeldas1D == 0 )
      {
        connec[ 6*fila + 1 ] = fila +aux6;
      }
      else
      {
        connec[ 6*fila + 1 ] = fila - 1;
      }
      if(fila == aux1 )
      {
        connec[ 6*fila + 2 ] = fila - aux7;
        aux1 = aux1+1;
        if(fila == (aux2-1) )
        {
            aux1 = aux2+ aux7;
            aux2 += NCeldas2D ;
        }
      }
      else
      {
        connec[ 6*fila + 2 ] = fila  + NCeldas1D ;
      }
      if( fila == aux3 )
      {
        connec[ 6*fila + 3 ] = fila + aux7;
        aux3 = aux3+1;
        if(fila == (aux4-1) )
            {
                aux3 = aux5;
                aux4 = aux3+NCeldas1D;
                aux5 += NCeldas2D;
            }
      }
      else
      {
        connec[ 6*fila + 3 ] = fila - NCeldas1D;
      }
      if(fila < (NCeldas1D-1)*NCeldas2D)
      {
        connec[ 6*fila + 4 ] = fila + NCeldas2D;
      }
      else
      {
        connec[ 6*fila + 4 ] = fila - (NCeldas1D-1)*NCeldas2D ;
      }
      if( fila < NCeldas2D)
      {
        connec[ 6*fila + 5 ] = fila + (NCeldas1D-1)*NCeldas2D ;
      }
      else
      {
        connec[ 6*fila + 5 ] = fila -NCeldas2D;
      }
    }
    
    //save connectivity matrix
    std::fill(connectivity.begin(),connectivity.end(),0);
    for(int j=0; j<dim_connec; j++)
    {
        connectivity[j] = connec[j];      
    }
    
    return;
}


// Fill head and list
// Be careful with indices in this algorithm, 0 has a special meaning in head and list,
// we do not identify the particles from 0 to N-1 but 1 to N. 
void make_linked_list3D(){
    
    std::fill(head.begin(), head.end(), 0);
    std::fill(list.begin(), list.end(), 0);

    int icell;   //Cell index in head
    double temppos[3]; //position of a particle

    for(int i=1; i<N+1; i++)
    {
        //change reference system
        temppos[0] = particles[3*(i-1)];
        temppos[1] = L-particles[3*(i-1)+1];
        temppos[2] = L-particles[3*(i-1)+2];
        
        //And reduce it to the primary box
        apply_pbc3D(temppos);
        //Compute the cell coordinates
        icell = getcell3D(temppos, i);
        make_configuration( i, icell);
        //Add particle to head and list
        list[i] = head[icell];
        head[icell] = i;
    }
    return;
} 

//Computes cell coordinates and cell index
int getcell3D(double *pos, int i){ 
    int icell; //Cell index for head
    int cell[3]; //mx, my, mz of a particle

    for(int j=0; j<3; j++)
    {
        cell[j] = (int) floor((*(pos+j)/l)+1.0);
    }
    icell = cell[0] + (cell[1]-1)*NCeldas1D+ (cell[2]-1)*NCeldas1D*NCeldas1D;
    
    //save
    molecules[3*i+0] = cell[0];
    molecules[3*i+1] = cell[1];
    molecules[3*i+2] = cell[2];
    
    return icell;
}


//create the configuration matrix 
void make_configuration(int i, int icell){
    int celda = (icell-1)*M + species[i];
    configuration[celda] +=1;
    return;
}
  
  
//initial rate matrix
void make_rateMatrix3D(){
    
    double K, q1, q2;            // M-1 is constant ratio K= k1/k2  
    
    for(int i=0; i<Nsv; i++)
    {
        
        if(k2==0)
        {
            q1 = 0;
            q2 = 0;
            K =0;
        }
        else{
            K = k1/k2;  
            q1 = change_rate3D(k1, D_V+D_A);   //association rate
            q2 = q1/K;                         //dissociation rate
        }
    
        double rate_r = 0, rate_s = 0;
    
        for(int es=0; es<M-2; es++)
        {
            rate_r +=  pvol*q1*configuration[M*i]*configuration[M*i+es+1] +  Molec2Molar*q2*configuration[M*i+es+2];
        }
        
        for(int sp=1; sp<M; sp++)
        {
            rate_s += (D_V/(l*l))*configuration[M*i+sp];
        }
        
        rate[3*i+0] = rate_r;
        rate[3*i+1] = (D_A/(l*l))*configuration[M*i] + rate_s;
        rate[3*i+2] = rate[3*i] + rate[3*i+1];
    }
   return;     
}     

//update rate matrix                
void update_rateMatrix3D(int s){
    
    double K, q1, q2;            // M-1 is constant ratio K= k1/k2 
    double rate_r = 0, rate_s = 0;
    
    if(k2==0)
    {
        q1 = 0;
        q2 = 0;
        K = 0;
    }
    else{
        K = k1/k2;  
        q1 = change_rate3D(k1, D_V+D_A);   //association rate
        q2 = q1/K;                         //dissociation rate
    }
    
    
    for(int es=0; es<M-2; es++)
    {
        rate_r +=pvol*q1*configuration[M*s]*configuration[M*s+es+1]+ Molec2Molar*q2*configuration[M*s+es+2];
    }
        
    for(int sp=1; sp<M; sp++)
    {
        rate_s += (D_V/(l*l))*configuration[M*s+sp];
    }
        
    rate[3*s+0] = rate_r;
    rate[3*s+1] = (D_A/(l*l))*configuration[M*s]+ rate_s;
    rate[3*s+2] = rate[3*s] + rate[3*s+1];
    
    return;
}     

//Create event queue array. First event time 
void make_event_q(){
    
    //Calculate t
    double t;     //t: time for the first reaction-diffusion event
    std::fill(event_queue.begin(), event_queue.end(), 0);
    std::fill(event.begin(), event.end(), 0);     
    
    for(int i=0; i<Nsv; i++)
    {
        if (rate[3*i+2]==0)
        {
            t = 1e5;
        }
        else{
            t = - log (Myrand())/ (rate[3*i+2]);
        }
        
        event[2*i+0] = i;
        event[2*i+1] = t;
        event_queue[2*i+0] = i;
        event_queue[2*i+1] = t;
    } 
    
    //event is sorted with increasing event time
    updateEvent();
    return;
} 

//Make Q-array. Subvolumes are ordering according increasing time
void make_q_matrix(){
    for(int ss=0; ss<Nsv; ss++)
    {
        q_array[ss]= event[2*ss]+1;
    }
    return;
}

//sort array of 2xn variables
void sort(int num_p){
    int pos, parent;                //heap positions
    vector<double> elem;            //for exchanging elements
    elem.resize(2);                //event queue array 
    std::fill(elem.begin(), elem.end(), 0); 
        
    pos = num_p;                           //insert at end
    parent = pos -1;                       //move up in heap
    
    //exchange parent/child
    while((pos > 0)&&((event[2*parent+1])>(event[2*pos+1])))
    {
        elem[1] = event[2*parent+1]; //time
        elem[0] = event[2*parent];   //ev
        event[2*parent+1] = event[2*pos+1];//time
        event[2*parent] = event[2*pos];   //ev
        event[2*pos+1] = elem[1]; //time
        event[2*pos] = elem[0];  //ev
        pos = parent;                //move up
        parent = pos-1;
    }
    return;
}


//reorder event matrix 
void updateEvent(){
    for(int s=0; s<Nsv; s++)
    {
        //for(int c=0; c<Nsv; c++){
        sort(s);
        //}
    }
    return;
}
    
    
// Reaction event
void reaction_event(int sv, int sw){
    
    //n_R: number of reaction
    int mol1, mol2, mol3, n_R;
    double K, q1, q2;            // M-1 is constant ratio K= k1/k2  
    
    if(k2==0)
    {
        q1 = 0;
        q2 = 0;
        K = 0;
    }
    else{
        K = k1/k2;  
        q1 = change_rate3D(k1, D_V+D_A);   //association rate
        q2 = q1/K;                         //dissociation rate
    }
    
    
    n_R = get_Reaction(sv, q1, q2);

    if(n_R <=(int(R/2)) )
    {
        if((configuration[M*sv]==0) || (configuration[M*sv+n_R]<0))
        {
	    printf("No reaction. No A \n");
            exit(EXIT_FAILURE);
        }
        else
        {
            r1++;
            r12 = r1+r2;
                        
            mol1 = get_molecule(sv, 0); 
            mol2 = get_molecule(sv, n_R); 
            mol3 = r1+N;
           
            update_moleculesR(sv,r12, mol1, mol2, mol3, 1, sw, n_R+1); 
            
            configuration[M*sv] = configuration[M*sv]-1;       
            configuration[M*sv+n_R] = configuration[M*sv+n_R]-1;   
            configuration[M*sv+n_R+1] = configuration[M*sv+n_R+1]+1; 
            
            //update rate matrix
            update_rateMatrix3D(sv); 
             
        }
    }
    else
    {
        
        if((configuration[M*sv+n_R-int(R/2)+1]==0) || (configuration[M*sv+n_R-int(R/2)+1]<0))
        {
            printf("No reaction. No C \n");
            exit(EXIT_FAILURE);
        }
        else
        {
            r2++;
            r12 =r1+r2;

            mol1 = 0;
            mol2 = 0;
            mol3 = get_molecule(sv, n_R-int(R/2)+1);
            update_moleculesR(sv,r12, mol1, mol2, mol3, 2, sw, n_R-int(R/2)); 
            
            //update configuration matrix
            configuration[M*sv] = configuration[M*sv]+1;       
            configuration[M*sv+int(n_R-R/2)] = configuration[M*sv+int(n_R-R/2)]+1; 
            configuration[M*sv+int(n_R-R/2+1)] = configuration[M*sv+int(n_R-R/2+1)]-1; 

            //update rate matrix
            update_rateMatrix3D(sv); 
        }
    }
        
    //Calculate next event (t_next)
    get_nextEvent(sv);
    return;
}


//diffusion event
void diffusion_event(int sv, int sw){
    
    double t_alpha = event_queue[2*sv+1];     // t_alpha 
    int mol, especie;
    double rs_old = 0.0, rs_o = 0.0;          // r_gamma+ s_gamma 
     
    especie = get_Diffusion(sv);
    
    int direccion;                               //next direction....0,1,2,3,4,5  array(6)
    int nsv=sv;                                  //nsv = new subvolume
    while (sv==nsv){
        direccion = (int) ((Myrand())*(2*dimension));
        nsv = connectivity[2*dimension*sv+direccion];
    }

    int celda = sv*M;      //old subvolume
    int celdanew = nsv*M;  //new subvolume
    
    
    if(configuration[celda+especie]==0)
    {
        printf("No particles:%d\n", especie+1);
        exit(EXIT_FAILURE);
    }
    else
    {
        //update configuration matrix 
        mol = get_molecule(sv, especie);
        configuration[celda+especie] -=1;
        configuration[celdanew+especie] +=1;
        
        //update rate matrix
        update_rateMatrix3D(sv);
        rs_old = rate[3*nsv+2];
        update_rateMatrix3D(nsv);
    }
    
   
    //move particles to the new volume
    move_particles(nsv, mol);
   
    //save data
    if(sw == 0){
        write_particles("timeParticlest.txt", event[1],mol);}
    if(sw == 1){
        write_particles("timeParticlest1.txt", event[1],mol);}
       
    //Calculate next reaction-diffusion event (t_next) 
    get_nextEventD(nsv, sv, rs_old, t_alpha);

    return;
}

//Get the particle that reacts 
int get_Reaction(int sv, double q1, double q2)
{
    std::vector<double> P_r; 
    P_r.resize(R, 0); 
    double result;
    int react;
    
   //Probability  
    for(int es=0; es<R/2; es++)
    {
        P_r[es]=   (pvol*q1*configuration[M*sv]*configuration[M*sv+es+1])/rate[3*sv];
        P_r[es+int(R/2)] =  (Molec2Molar*q2*configuration[M*sv+es+2])/rate[3*sv];
    }
    
    auto res = std::max_element(std::begin(P_r), std::end(P_r));
    
    if (std::end(P_r)!=res){
        result= *res;}
        
    for(int s=0;s<R;s++){ 
        if(P_r[s] ==*res){react = s+1;}
    }
    
    return react;
}


//Get the particle that diffuses
int get_Diffusion(int sv)
{
    std::vector<double> P_d; 
    P_d.resize(M, 0); 
    double result;
    int diff;
    
    //Probability to diffuse  
    for(int es=1; es<M; es++)
    {
        // Probability to diffuse A
        P_d[0]= (D_A/(l*l))*configuration[M*sv]/rate[3*sv+1]; 
        // Probability to diffuse V, VA...
        P_d[es] = (D_V/(l*l))*configuration[M*sv+es]/rate[3*sv+1];
    }
    
    auto res = std::max_element(std::begin(P_d), std::end(P_d));
    
    if (std::end(P_d)!=res){
        result= *res;}
        
    for(int s=0;s<M;s++){ 
        if(P_d[s] ==*res){diff = s;}
    }
    
    return diff;
}


//Calculate next reaction-diffusion event (t_next)  and Update event matrix
void get_nextEventD(int nsv, int sv, double rs_old, double t_alpha){
    
    double tnew;
    double tnewn;
    if(rate[3*sv+2]==0)
    {
        tnew = INFINITY;
    }
    else
    {
        tnew = -log(Myrand())/rate[3*sv+2];
    }
    
    //time for event matrix (subvolume gamma)
    if(rate[3*nsv+2]==0)
    {
        tnewn = INFINITY;
    }
    else
    {
        if(rs_old == 0){
            tnewn = -log(Myrand())/rate[3*nsv+2];
        }
        else{
            //tnewn = (rs_old/rate[3*nsv+2])*(event_queue[2*nsv+1]-t_alpha); //time of Gibson-Bruck
            tnewn = -log(Myrand())/rate[3*nsv+2];
        } 
    }

    
    //Update event matrix
    event_queue[2*sv+1] = tnew + t_alpha;        //Update event matrix (subvolume=lambda) 
    event_queue[2*nsv+1] = tnewn + t_alpha;     //Update event matrix (subvolume=gamma). Gibson-Bruck  
    
    event[1] += tnew;           
    for(int m=0; m<Nsv; m++)
    {
        if((nsv)==event[2*m])
        {
            event[2*m+1] = event_queue[2*nsv+1]; 
            break;
        }
    }
    
    //reorder event matrix 
    updateEvent();
    
    return;
}

//Calculate next reaction-diffusion event (t_next)  and Update event matrix
void get_nextEvent(int sv){
    
    double tnew;
    if(rate[3*sv+2]==0)
    {
        tnew = 1e5;
    }
    else
    {
        tnew = - log(Myrand())/rate[3*sv+2];
    }
    
    //Update event matrices
    event_queue[2*sv+1] = event_queue[2*sv+1] + tnew;      
    event[1] = event[1] + tnew;
    
    updateEvent();
    
    return;
}


//get the molecule that moves
int get_molecule(int sv, int especie){
    
    //choose random particle to move
    int Nmol = configuration[M*sv+especie];
    int j = (int) (Myrand()*Nmol);                          
    int i, p, p0, p1, p2, p3;
    vector<int> mol_sv;              
    mol_sv.resize(Nmol,0);             
    p = 0;
    p0 = head[sv+1];
    
        
    //get molecule from head
    if(species[p0]==especie){
        mol_sv[p] = p0;
        p = p+1;
    }
    
    //get molecule from list
    if(head[sv+1]!= 0)
    {
        p1 = list[head[sv+1]];
        if(species[p1]==especie)
        {
            mol_sv[p] = p1;
            p = p+1;
        }
       
        while(list[p1]!= 0)
        {
            p2 = list [p1];
            if(species[p2]==especie){
                mol_sv[p] = p2;
                p = p+1;
            }
            p1 = p2;
        }
    }
    
    i = mol_sv [j];
    
    return i; 
}

//save position of molecules  
void move_particles(int nsv, int i){
    
    //save molecule  
    for(int jj=0; jj<3; jj++)
    {
        molecules[3*i+jj] = positionCell[(dimension+1)*nsv+jj+1];
    }

    //update head and list
    update_headList();

    return;
} 


//Update head array and list array
void update_headList(){   
    
    std::fill(head.begin(), head.end(), 0);
    std::fill(list.begin(), list.end(), 0);

    int icell;   
    int temppos[3]; 
    for(int i=1; i<N+N_C+1; i++)
    {
        temppos[0] = molecules[3*i];
        temppos[1] = molecules[3*i+1];
        temppos[2] = molecules[3*i+2];
       
        if(species[i]!=M+1)
        {
            //Compute the cell coordinates
            icell = temppos[0] + (temppos[1]-1)*NCeldas1D+ (temppos[2]-1)*NCeldas1D*NCeldas1D;
            list[i] = head[icell];
            head[icell] = i;
        }
    }
   
    return;
}


//determine positions of subvolumes
void make_positionCell(){

    int a=0;
    int b=0;
    int c=0;
    for(int i=0; i<Nsv; i++)
    {        
        positionCell[(dimension+1)*i] = i;
        if(a==NCeldas1D){
            a = 1;
            b = b+1;
            positionCell[(dimension+1)*i+1] = a;
            }
        else{
            positionCell[(dimension+1)*i+1] = a+1;
            a = a+1;
            }
        if(b==NCeldas1D){
            b=0;
            positionCell[(dimension+1)*i+2]= b+1;
            }
        else{
            positionCell[(dimension+1)*i+2]= b+1;
            }
        if(i!=0 && (i+1)%NCeldas2D==0.0){
            c=c+1;
            positionCell[(dimension+1)*i+3]= c;
            }
        else{
            positionCell[(dimension+1)*i+3]= c+1;
            }
    } 
    
    return;
}

//save molecules of reactions and update lists 
void update_moleculesR(int sv, int r, int mol1, int mol2, int mol3, int save , int sw, int esp){

    //save molecules
    if(save ==1){

        for(int j=0; j<3; j++)
        {
            //save positions x y z de moleculesC
            moleculesC[6*r+j] = positionCell[(dimension+1)*sv+j+1];
            molecules[3*mol1+j] = positionCell[(dimension+1)*sv+j+1];
            molecules[3*mol2+j] = positionCell[(dimension+1)*sv+j+1];
            molecules[3*mol3+j] = positionCell[(dimension+1)*sv+j+1];
        }
        //Update molecules from reactions
        moleculesC[6*r+3] = mol1;
        moleculesC[6*r+4] = mol2;
        moleculesC[6*r+5] = mol3;     
        //Update species
        species[mol1] = M+1;                 
        species[mol2] = M+1;                 
        species[mol3] = esp;                 
    }
    
    if(save == 2)
    {
        
        for(int m=1; m<N_C+1; m++)
        {   
            if(moleculesC[6*m+5]==mol3){
                mol1 = moleculesC[6*m+3];
                mol2 = moleculesC[6*m+4];
                break;
            }
        } 
        
   
        for(int j=0; j<3; j++)
        {
            moleculesC[6*r+j] = positionCell[(dimension+1)*sv+j+1];
            
            molecules[3*mol1+j] = positionCell[(dimension+1)*sv+j+1];
            molecules[3*mol2+j] = positionCell[(dimension+1)*sv+j+1];
            molecules[3*mol3+j] = positionCell[(dimension+1)*sv+j+1];
        }
        
        //update molecules from reactions
        moleculesC[6*r+3] = mol1;
        moleculesC[6*r+4] = mol2;
        moleculesC[6*r+5] = mol3;
        
        //Update species
        species[mol1] = 0;                   
        species[mol2] = esp;                 
        species[mol3] = M+1;                 
        
    }

    //update head and list
    update_headList();
    
    if(sw == 0){
        if(save==1){
            write_particles("timeParticlest.txt", event[1],mol3);  
        }
        
        if(save==2){
            write_particles("timeParticlest.txt", event[1],mol1);
            write_particles("timeParticlest.txt", event[1],mol2);
        }
    }
    
    if(sw == 1){
        if(save==1){
            write_particles("timeParticlest1.txt", event[1],mol3);}
        if(save==2){
            write_particles("timeParticlest1.txt", event[1],mol1);
            write_particles("timeParticlest1.txt", event[1],mol2);}
    }
    return;
}


//total number of particles of each specie:  t-A-V-A1V-A2V-...
void write_particlestimeEvolution(const char* name, double t){    
    static ofstream out1(name);
    vector<double> part;
    part.resize(M+1,0);
    
    part[0] = t;
    
    for(int ss=0; ss<Nsv; ss++)
    {
        for(int p=0; p<M; p++){
            part[p+1] += configuration[M*ss+p];
        }
    }
    
    for(int i =0; i<M+1; i++)
    {
        out1 << part[i] << "  "  ;
    }
    out1 <<"\n";

    out1.flush();
  
}

//All particles at time t  : "particle-t-x-y-z-especie"
void write_allparticles(const char* filename, double t){   
    static ofstream gout(filename);
    int positions[4];
    for(int i=1; i<N+r12+1; i++)
    {
        for(int j=0; j<dimension+1; j++)
        {
        positions[j] = molecules[dimension*i+j];
        positions[3] = species[i];
    }
    gout << i << " " << t << " " << positions[0] << " " << positions[1] << " " << positions[2] << " " << positions[3] << " " <<"\n";
   }
   gout.flush();
}

//the particle in an event at time t  : "particle-t-x-y-z-especie"
void write_particles(const char* name, double t, int i){   
    static ofstream gout(name);
    int positions[4];
    for(int j=0; j<dimension+1; j++)
    {
        positions[j] = molecules[dimension*i+j];
        positions[3] = species[i];
    }
    gout << i << " " <<fixed<<setprecision(10)<< t << " " << positions[0] << " " << positions[1] << " " << positions[2] << " " << positions[3] << " " <<"\n";
    gout.flush();
}

//Periodic Boundary Conditions
void apply_pbc3D(double *r){
    for(int i=0; i<3; i++)
    {
        while( r[i] > L )
        {
            r[i] -= L;
        }
        while( r[i] < 0 )
        {
            r[i] += L;
        }
    }
    return;
}


//change reaction rate constant (mesoscopic limit) 
double change_rate3D(double ka, double D){
    
    //Degree of diffusion control 
    double alpha;
    double pi = 3.14159265358979323846;
    double k, ka2;
    ka2 = ka*Molar2Molec;
    k = (4*pi*rho*D*ka2)/ (4*pi*rho*D - ka2);                     // k microscopic
    alpha = k / (4*pi*rho*D);   // degree of diffusion control 3D
    
    //spatial discretization 
    double h = l*(cbrt(21/(4*pi)))-rho;
    double R = rho+h;
    //unit-less measure for spatial discretization
    double beta = rho/(rho+h);
    double q;
    
    //mesoscopic rate constant
    q = k/ (1+alpha*(1-beta)*(1-0.58*beta));
    return q;
}


//write event matrix
void writeEvent(){   
	FILE *fp;
	fp = fopen("event.txt", "w+");
	
	for(int j=0; j<Nsv; j++)
    {
		fprintf(fp, "%g %g \n", event_queue[2*j], event_queue[2*j+1]);
    }
    fclose(fp);
	
	return;
}


//write molecules
void writeMolecules(){   
    
    FILE *fp;
    fp = fopen("molecules.txt", "w+");
    
    for(int j=0; j<N+r1+1; j++)
    {
		fprintf(fp, "%d %d %d \n", molecules[3*j], molecules[3*j+1], molecules[3*j+2]);
    }
    fclose(fp);
    return;
}

//write species
void writeSpecies(){   

    FILE *fp;
    fp = fopen("species.txt", "w+");
    
    for(int j=0; j<N+r1+1; j++)
    {
		fprintf(fp, "%d\n",species[j]);
    }
	
    fclose(fp);
    return;
}


//write configuration
void writeConfiguration(){   
    
    FILE *fp;
    fp = fopen("configuration.txt", "w+");
    
    for(int j=0; j<Nsv; j++)
    {
        for(int tipo_particula = 0; tipo_particula < M; tipo_particula++)
        {
		  fprintf(fp, "%d ", configuration[M*j + tipo_particula ]);
        }
        fprintf(fp, "\n");
    }
	
    fclose(fp);
    return;
}

//read event
void readEvent(){
    
    FILE *fp;
	double values[2*Nsv];
    unsigned int i;  
    std::fill(event_queue.begin(), event_queue.end(), 0);
    std::fill(event.begin(), event.end(), 0);  

fp = fopen("event.txt", "r"); 

    if(fp == NULL){

        perror("Error");
        exit(EXIT_FAILURE);
    }
    else
    {          
        for(i = 0; i < 2*Nsv; ++i){
        int m = fscanf(fp, "%lf",&values[i]);
        }
    
        for(int j = 0; j < 2*Nsv; ++j){
                event[j] = values[j];
                event_queue[j] = event[j];
        }
    
    }
    
    fclose(fp);

    return;
}

//read Molecules
void readMolecules(){
    
    FILE *fp;
	double values[3*(N+r1+1)];
    unsigned int i; 
    
	fp = fopen("molecules.txt", "r"); 

    if(fp == NULL){
        perror("Error");
        exit(EXIT_FAILURE);
    }
    else
    {           
        for(i = 0; i < 3*(N+r1+1); ++i){
        int m = fscanf(fp, "%lf",&values[i]);
        }
    
        for(int j = 0; j < 3*(N+r1+1); ++j){
            molecules[j]= (int) values[j];        
        }
    
    }
    fclose(fp);
    
    return;
}

//read species
void readSpecies(){
    
    FILE *fp;
	double values[N+r1+1];
    unsigned int i; 
    int na=0, nv=0, nav=0;
    
	fp = fopen("species.txt", "r"); 
    
    if(fp == NULL){
        perror("Error");
        exit(EXIT_FAILURE);
    }
    else
    {          
        for(i = 0; i < N+r1+1; ++i){
        int m = fscanf(fp, "%lf",&values[i]);
        }
    
        for(int j = 0; j < N+r1+1; ++j){
            species[j]= (int) values[j]; 
        }
    }
    
    fclose(fp);

    return;
}


//read configuration
void readConfiguration(){   
    
    FILE *fp;
	double values[M*Nsv];
    unsigned int i;   
    
	fp = fopen("configuration.txt", "r");

    if(fp == NULL){
        perror("Error");
        exit(EXIT_FAILURE);
    }
    else
    {          
        for(i = 0; i < M*Nsv; ++i){
        int m = fscanf(fp, "%lf",&values[i]);
        }
    
        for(int j = 0; j < M*Nsv; ++j){
            configuration[j]=  (int) values[j];        
        }
    
    }
    
    fclose(fp);

    return;
}

void writeInfo(){   

    FILE *fp;
    fp = fopen("info.txt", "w+");
    
	fprintf(fp, "%d \n",r1);
    fprintf(fp, "%d \n",r2);
    fprintf(fp, "%d \n",r12);
    fprintf(fp, "%f \n",tclock);
	
    fclose(fp);
    return;
}


//read Info
void readInfo(){   
    
    FILE *fp;
	double values[4];  
    unsigned int i;
    
	fp = fopen("info.txt", "r");

    if(fp == NULL){
        perror("Error");
        exit(EXIT_FAILURE);
    }
    else
    {   
        for(i = 0; i < 3; ++i){
            int m = fscanf(fp, "%lf",&values[i]);}
            
        int react =  (int) values[0];
        r1 = react;
        r2 =  (int) values[1];
        r12 =  (int) values[2];
    }
    fclose(fp);

    return;
}

void writeMolC(){   

    FILE *fp;
    fp = fopen("molC.txt", "w+");
	
    for(int j=0; j<r12+1; j++)
    {
		fprintf(fp, "%d %d %d %d %d %d \n", moleculesC[6*j], moleculesC[6*j+1], moleculesC[6*j+2], moleculesC[6*j+3], moleculesC[6*j+4], moleculesC[6*j+5]);
    }
	
    fclose(fp);
    return;
}


void readMolC(){   
    
    FILE *fp;
	double values[6*(r12+1)];
    unsigned int i; 
    
	fp = fopen("molC.txt", "r");

    if(fp == NULL){

        perror("Error");
        exit(EXIT_FAILURE);
    }
    else
    {   
        for(i = 0; i < (6*(r12+1)); ++i){
        int m = fscanf(fp, "%lf",&values[i]);
        }
    
        for(int j = 0; j < (6*(r12+1)); ++j){
            moleculesC[j]= (int) values[j]; 
        }
    }
    
    fclose(fp);

    return;
}

//write positions of C "t-x-y-z-V-A-VA-r1"
void write_timeEvolutionC(const char* name , double t){   
    static ofstream gout(name);
    double positions[6];
    for(int j=0; j<dimension+3; j++)
    {
            positions[j] = moleculesC[6*r12+ j];
    }
    gout << t << " " << positions[0] << " " << positions[1] <<  " " << positions[2] << " " << positions[3] << " " << positions[4] <<  " " << positions[5] << " " << r1 << " " <<"\n";
    
    gout.flush();
}

//read initial parameters
void readParameters(){   
    
    FILE *fp;
	double values[13];
    unsigned int i;   
    
	fp = fopen("parameters.txt", "r");

    if(fp == NULL){
        perror("Error");
        exit(EXIT_FAILURE);
    }
    else
    {          
        for(i = 0; i < 13; ++i){
        double m = fscanf(fp, "%lf",&values[i]);
        }

        for(int j = 0; j < 13; ++j){
            parameters[j]=  (double) values[j];        
        }   
    }
    
    fclose(fp);
    return;
}

//Write a summary of simulation
void writeSummary(){   

    FILE *fp;
    fp = fopen("summarySimulation.txt", "w+");
    
    fprintf(fp, "No particles of type A: %d\n", N_A);
    fprintf(fp, "No particles of type V: %d\n", N_V);
    fprintf(fp, "No of reactions  %d \n",R/2);
    fprintf(fp, "kon (s-1M-1): %0.2f \n",k1);
    fprintf(fp, "koff (s-1):%0.2f \n",k2);

    fclose(fp);
    return;
}

//random number between 0 and 1
double Myrand(){
return ((double) ran2(seed));    
}

