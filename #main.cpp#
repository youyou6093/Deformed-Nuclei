#include <chrono>
#include "vector_compute.h"
#include "dirac_oscillator.h"
#include <iostream>
#include "generate_state.h"
#include <unordered_map>
#include <string>
#include <fstream>
#include "integrator.h"
#include "generate_matrix.h"
#include "simps.h"
#include "symm.h"
#include "Solution.h"
#include "Density.h"
#include "Green-method.h"
#include "effective_density.h"
#include "utility.h"
#include <algorithm>
using namespace std;
#define MAX_KAPPA 7

/*--------------------set up global varialbes-----------------------------*/

/*constants*/
const double PI=atan(1.0)*4;
string parameters = "parameter.txt";  // the file contains the model parameters
double my_tolerance=1e-3;   //tolerance of the integrator

/*model parameters, default values are FSUGOLD*/
double hbarc=197.326;
double ms=491.500/hbarc;
double mv=782.5/hbarc;
double mp=763.0/hbarc;
double mg=0.00000001/hbarc;
double gs=112.199551;
double gv=204.546943;
double gp=138.470113;
double gg=4*PI/137.0;
double lambdas=0.0;
double lambdav=0.030;
double lambda=0.023762;
double ks=0.0600;
double ka=1.420333/hbarc;


/* Computation constants */
double b=2.4,rmin=0.0,rmax=20.0;
int N=801;
/*Max value of L CHANNEL, it is initialized to 1( nothing at all! ).
 However, it will change after 50(or other values) iterations */
int max_L = 4;
/*The potential I want to add, initialized to 0 but
 will change to 1(or other value) after some iterations */
int Potential_channel = 1;
/*After a few iterations it will also change */
double Deformation_parameter = 0.05;
double h=(rmax-rmin)/(N-1);       //grid size
int proton_number, neutron_number;
vector<double> start ={450.0,400.0,5.0,0.0}; //initial value of Phi W B A in the wood-saxon form
vector<double> fx;                //X value array
vector<double> empty;             //empty vector, will be initialized later
unordered_map<string, vector<vector<double>> > States_map;   //store the infos about the wavefunctions when given states
unordered_map<string, double> Energy_map;                    //store the energies for every state
unordered_map<string, double > Angular_map;                  //Stroe the value A in my notes
/*---------------------------------finishing global varialbes-------------------------------------*/


#include "preprocess.h"









int main(){
    /*test the time*/
    chrono::steady_clock::time_point tp1 = chrono::steady_clock::now();
    unordered_map<string, double>::iterator Angular_ptr;
    unordered_map<string, vector<vector<double>>>::iterator States_ptr;
    vector<vector<double>> Phi,W,B,A,dens,denv,denp,den3,EFF_Phi,EFF_W,EFF_B,EFF_A;
    vector<vector<double>> scalar_p,vector_p,scalar_n,vector_n;
    vector<double> Potential,flat; //flat is the flatted matrix
    vector<State> States_m;        //the basis quantum number of a specific m
    matrix_diag diag = matrix_diag();      //diagnolization class
    vector<eig2> occp,occn;        //occupied states of protons and neutrons
    vector<eig2> occp_raw,occn_raw;      //raw solution of matrix
    vector<eig2> temp_solution;
    /*determine the Z and N;*/
    cout << "Z=: ";
    cin  >> proton_number;
    cout << "N=: ";
    cin >> neutron_number;
    cout << "Deformation_parameter=: ";
    cin >> Deformation_parameter;
    preprocessing();                    //create the states_map and angular_map
    preprocessing_2(Phi, W, B, A);      //initialize the PHI W B A potential
    /*Initialize the densities to 0*/
    for(int i=0;i<max_L;i++){
        dens.push_back(empty);
        denv.push_back(empty);
        denp.push_back(empty);
        den3.push_back(empty);
    }
    /*Initially, there is no dipole potential added*/
    for(int i=0;i<N;i++)               
        Potential.push_back(0.00); 
    EFF_Phi=dens; EFF_W=dens;EFF_A=dens;EFF_B=dens;   //set everything to zero;
    /*determine the range of m*/
    double min_m = -magic(max(proton_number,neutron_number));
    double max_m = -min_m;
    
    cout<<"finishing processing"<<endl;
    cout << "Deformation_parameter = " << Deformation_parameter << endl;
    /*start iteration */
    vector<Solution> Final_occp,Final_occn;                 //Final solution for 1 iteration
    
    
    for(int ite=0;ite<50;ite++){
        chrono::steady_clock::time_point tpold = chrono::steady_clock::now();
        occp_raw.clear();
        occn_raw.clear();
        occp.clear();
        occn.clear();            //make sure they are empty
        for(int i=0;i<N;i++)                //initialize the excitation term
            Potential[i] = fx[i] * Deformation_parameter;
        /*because every iteration I update the meson potentials so I need to recompute the scalar
         and vector potentials*/
        generate_potential(Phi, W, B, A, Potential, scalar_n, vector_n, Potential_channel, 0);   //neutron
        generate_potential(Phi, W, B, A, Potential, scalar_p, vector_p, Potential_channel, 1);   //proton
        
        /* solve the dirac equation for one iteration */
        for(double m = min_m ; m < max_m + 1 ; m++){
            States_m=generate_statesm(m, max_L);        //states for this m
            //generate_full_matrix(scalar_n, vector_n, States_m, m);
            vector<vector<double>> M = generate_full_matrix(scalar_n, vector_n, States_m, m);
            if (M.size()!=States_m.size()) cout<<"error!!"<<endl;
            /*another test */
            flat=flat_matrix(M);
            diag=matrix_diag(flat,int(States_m.size()));    //diagnolize the matrix, add all the eigenvalues and eigenvectors into M_matrix
            diag.get_results();                             //diag.results is a matrix full of eigvalues and eigenvectors.
            temp_solution=get_temp_solution(diag.results, m);   //add m to the results
            occn_raw.insert(occn_raw.end(), temp_solution.begin(), temp_solution.end());
            M=generate_full_matrix(scalar_p, vector_p, States_m, m);    //get the matrix for the specific m,proton
            flat=flat_matrix(M);
            diag=matrix_diag(flat,int(States_m.size()));
            diag.get_results();
            temp_solution=get_temp_solution(diag.results, m);
            occp_raw.insert(occp_raw.end(),temp_solution.begin(),temp_solution.end());
        }
        
        get_solution(occp_raw, occn_raw, occp, occn);            //get sorted occ state
        Final_occn=get_solutions_object(occn);                   //get solution object
        Final_occp=get_solutions_object(occp);
    
        generate_density(occn,occp,dens,denv,denp,den3); //calculate all the density based on the occupied states
        
        for(int i = 0; i < max_L; i++)
            cout << "channel=" << i << ':' << dens[i][0] << ' ' << denv[i][0] << ' ' << den3[i][0] << ' ' << denp[i][0] << endl;
        /* compute the effective density
         solve the klein-gordon equation */
        update_potential(EFF_Phi, EFF_B, EFF_A, EFF_W, Phi, W, B, A, dens, denv, den3, denp);
        for(int i = 0; i < max_L; i++)
            cout <<"channel="<<i<<':'<< Phi[i][0] <<' '<< W[i][0]<<' '<<B[i][0]<<' '<<A[i][0]<<endl;
        /*check how many occ states we found for each iteration*/
        cout<<occp.size()<<' '<<occn.size()<<endl;
        /*get energy*/
        cout<<"E/A="<<compute_energy(occp, occn, Phi, W, B, A, dens, denv, den3, denp)<<endl;
        
        chrono::steady_clock::time_point tpnew = chrono::steady_clock::now();
        chrono::steady_clock::duration duration_in_ites = tpnew - tpold;
        cout << "Time_used_in_iteration " <<  ite << " = " <<  chrono::duration_cast<chrono::seconds>(duration_in_ites).count() << endl;
    } //end the iterations
    
    /* option part ,
     out put data */
    for(int i = 0 ;i < max_L; i++){
        ofstream outfile;
        ofstream outfile2;
        outfile.open("density" + to_string(i) + ".txt");
        outfile2.open("potential" + to_string(i) + ".txt");
        for(int j = 0; j < N; j++){
            outfile << fx[j] << ' ' << dens[i][j] << ' ' << denv[i][j] << ' ' << den3[i][j] << ' ' << denp[i][j] << endl;
            outfile2 << fx[j] << ' ' << Phi[i][j] << ' ' << W[i][j] << ' ' << B[i][j] << ' ' << A[i][j] << endl;
        }
        outfile.close();
        outfile2.close();
    }
    
    /* get time */
    chrono::steady_clock::time_point tp2 = chrono::steady_clock::now();
    chrono::steady_clock::duration d = tp2-tp1;
    cout <<"Total time used:"<<chrono::duration_cast<chrono::seconds>(d).count()<<endl;
    cout <<"maxL = "<<max_L<< endl;
    cout <<"Deformation parameter = "<< Deformation_parameter << endl;
}









/* the newly implemented part, when
 The iteration is bigger thant some presetted time,
 I will add dipole potential and add more L channel*/
//        if (ite > -1){
//            max_L = 3;
//            Deformation_parameter =  0.1 ;
//            Potential_channel = 1;
//            for (int i = 1 ; i < max_L; i++){
//                Phi.push_back(empty);
//                W.push_back(empty);
//                B.push_back(empty);
//                A.push_back(empty);
//                dens.push_back(empty);
//                denv.push_back(empty);
//                denp.push_back(empty);
//                den3.push_back(empty);
//            }
//            EFF_Phi=dens; EFF_W=dens;EFF_A=dens;EFF_B=dens;
//        }









