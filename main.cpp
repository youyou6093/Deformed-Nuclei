/*
 The added potential right now is r Y10 and in units Mev.
*/

//#include <cstdint>
#include <chrono>
#include "vector_compute.h"
#include "generate_state.h"
#include "dirac_oscillator.h"
#include <iostream>
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
#define ITENUM 200

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
/*Max value of L CHANNEL*/
int max_L = 1;
/*The potential I want to add*/
int Potential_channel = 1;
double Deformation_parameter = 0.00;
double h=(rmax-rmin)/(N-1);       //grid size
int proton_number, neutron_number;
vector<double> start ={450.0,400.0,5.0,0.0}; //initial value of Phi W B A in the wood-saxon form, units in Mev
vector<double> fx;                //x axis
vector<double> empty;             //empty vector, will be initialized later
unordered_map<string, vector<vector<double>> > States_map;   //store the infos about the wavefunctions when given states
unordered_map<string, double> Energy_map;                    //store the energies for every state
// unordered_map<string, double > Angular_map;                  //Stroe the value A in my notes
int max_k; //max kappa when I want to construct the basis 
/*---------------------------------finishing global varialbes-------------------------------------*/
#include "preprocess.h"







/*body of the function*/
int main(int argc, char ** argv){
    
    
    /*test part*/
    
    cout << Angular_depedencek(14,14,27,0) << endl;
//    cout << Angular_depedencek(14,14,-27,0) << endl;
    
    
    
    
    
    
    
    
    chrono::steady_clock::time_point tp1 = chrono::steady_clock::now();  //current time
    unordered_map<string, double>::iterator Angular_ptr;
    unordered_map<string, vector<vector<double>>>::iterator States_ptr;
    vector<vector<double>> Phi,W,B,A,dens,denv,denp,den3,EFF_Phi,EFF_W,EFF_B,EFF_A;
    vector<vector<double>> dens_old, denv_old, denp_old, den3_old;
    vector<vector<double>> scalar_p,vector_p,scalar_n,vector_n;
    vector<double> Potential,flat; //flat is the flatted matrix
    vector<State> States_m;        //the basis quantum number of a specific m
    matrix_diag diag = matrix_diag();      //diagnolization class
    vector<eig2> occp,occn;        //occupied states of protons and neutrons
    vector<eig2> occp_raw,occn_raw;      //raw solution of matrix
    vector<eig2> temp_solution;

    //Command line part
    proton_number = atoi(argv[1]);
    neutron_number = atoi(argv[2]);
    Deformation_parameter = stod(argv[3]);
    max_L = atoi(argv[4]);
    max_k = atoi(argv[5]);

    preprocessing();                    //Prepare the states Hash Map
    /*initialize the PHI W B A potential
      the potentials here are all in Mev*/
    preprocessing_2(Phi, W, B, A);
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
    /*Set the effecitive density in KG part to be the same as
      original density*/
    EFF_Phi=dens; EFF_W=dens;EFF_A=dens;EFF_B=dens;
    /*Determine the range of m, just roughly determine*/
    int min_m = -magic(max(proton_number,neutron_number));
    int max_m = -min_m;
    cout<<"Finishing processing"<<endl;
    cout << "Z = " << proton_number << " " << "N = " << neutron_number << endl;
    cout << "Deformation_parameter = " << Deformation_parameter << endl;
    /*The occupied states from solving */
    vector<Solution> Final_occp,Final_occn;                 //Final solution for 1 iteration
    
    
    for(int ite=0;ite<ITENUM;ite++){
        chrono::steady_clock::time_point tpold = chrono::steady_clock::now();
        /*Make sure OCCs are empty*/
        occp_raw.clear();
        occn_raw.clear();
        occp.clear();
        occn.clear();
        /* The potential term added to my model*/
        // for(int i=0;i<N;i++){
        // 	if (fx[i] <= 6)
        //     	Potential[i] = fx[i] * Deformation_parameter * sqrt(3.0 / 4 / PI);
        //     else
        //     	Potential[i] = 0.0;
        // }
        for(int i=0;i<N;i++)
        	Potential[i] = fx[i] * Deformation_parameter * sqrt(3.0 / 4 / PI);

        // * exp(- 0.4 * fx[i])
        /*Every iteration it updates the meson potentials
          ,so it need to recompute the scalar
           and vector potentials*/
        
        /*those two functions are in preprocess.h*/
        generate_potential(Phi, W, B, A, Potential, scalar_n, vector_n, Potential_channel, 0);   //neutron   preprocess.h
        generate_potential(Phi, W, B, A, Potential, scalar_p, vector_p, Potential_channel, 1);   //proton
        
        /* Solve the dirac equation*/
        for(int m = min_m ; m < max_m + 1 ; m+=2){      //The m here is actually 2 m
            States_m=generate_statesm(m, max_L);        //States for this m
//            cout << "States number is : " << States_m.size() << endl;
            vector<vector<double>> M = generate_full_matrix(scalar_n, vector_n, States_m, m);
            if (M.size()!=States_m.size()) cout<<"Error!!"<<endl;
            flat=flat_matrix(M);
            /*diagnolize the matrix, add all the eigenvalues and eigenvectors into M_matrix*/
            diag=matrix_diag(flat,int(States_m.size()));
            diag.get_results();
            /*diag.results is a matrix full of eigvalues and eigenvectors.
              I need to add m to the result to form the basis*/
            temp_solution = get_temp_solution(diag.results, m);
            /*Insert all solutions for this m*/
            occn_raw.insert(occn_raw.end(), temp_solution.begin(), temp_solution.end());
            /*get the matrix for the specific m,proton*/
            M=generate_full_matrix(scalar_p, vector_p, States_m, m);
            flat=flat_matrix(M);
            diag=matrix_diag(flat,int(States_m.size()));
            diag.get_results();
            temp_solution = get_temp_solution(diag.results, m);
            occp_raw.insert(occp_raw.end(),temp_solution.begin(),temp_solution.end());
        }
        
        /* Form the solution for occupied state*/
        get_solution(occp_raw, occn_raw, occp, occn);            //get sorted occ state, utility.h
        Final_occn=get_solutions_object(occn);                   //get solution object, density.h
        Final_occp=get_solutions_object(occp);
        /* Form the densities */
        generate_density(occn,occp,dens,denv,denp,den3); //calculate all the density based on the occupied states, density.h
        
//        if(ite > 0){
//            for(int i = 0; i < max_L; i++){
//                 for(int j = 0; j < N; j++){
//                     dens[i][j] = 3./4 * dens_old[i][j] + 1.0/4 * dens[i][j];
//                     denv[i][j] = 3./4 * denv_old[i][j] + 1.0/4 * denv[i][j];
//                     den3[i][j] = 3./4 * den3_old[i][j] + 1.0/4 * den3[i][j];
//                     denp[i][j] = 3./4 * denp_old[i][j] + 1.0/4 * denp[i][j];
//                 }
//             }
//         }
        
        dens_old = dens;
        denv_old = denv;
        den3_old = den3;
        denp_old = denp;
        
        for(int i = 0; i < max_L; i++)
            cout << "Channel=" << i << ':' << dens[i][0] << ' ' << denv[i][0] << ' ' << den3[i][0] << ' ' << denp[i][0] << endl;
        /* compute the effective density
         solve the klein-gordon equation
         update the potentials*/
        update_potential(EFF_Phi, EFF_B, EFF_A, EFF_W, Phi, W, B, A, dens, denv, den3, denp);
        for(int i = 0; i < max_L; i++)
            cout <<"channel="<<i<<':'<< Phi[i][50] <<' '<< W[i][50]<<' '<<B[i][50]<<' '<<A[i][50]<<endl;
        /*check how many occ states we found for each iteration*/
        cout<<occp.size()<<' '<<occn.size()<<endl;
        
        /* output the energy*/
        // for(int ii = 0; ii < occp.size(); ii ++){
        //     eig2 test = occp[ii];
        //     Solution temp = Solution(test);
        //     cout << ii << " energy for occp = " << temp.energy << ' ' << "m = " << temp.m << endl;
        //     temp.print_eigenvectors();
        // }
        // for(int ii = 0; ii < occn.size(); ii ++){
        //     eig2 test = occn[ii];
        //     Solution temp = Solution(test);
        //     cout << ii << " energy for occn = " << temp.energy << ' ' << "m = " << temp.m << endl;
        //     temp.print_eigenvectors();

        // }
        
        
        for(int i = 0 ;i < max_L; i++){
            if((ite % 5) == 0){

            ofstream outfile;
            ofstream outfile2;
            outfile.open("output/density" + to_string(i) + ' ' + to_string(ite) + ".txt");
            // outfile2.open("output/potential" + to_string(i) + ' ' + to_string(ite) +  ".txt");
            for(int j = 0; j < N; j++){
                outfile << fx[j] << ' ' << dens[i][j] << ' ' << denv[i][j] << ' ' << den3[i][j] << ' ' << denp[i][j] << endl;
                // outfile2 << fx[j] << ' ' << Phi[i][j] << ' ' << W[i][j] << ' ' << B[i][j] << ' ' << A[i][j] << endl;
            }
            outfile.close();
            outfile2.close();
        }
        }
        
        
        /*get energy*/
        cout<<"E/A="<<compute_energy(occp, occn, Phi, W, B, A, dens, denv, den3, denp)<<endl;
        chrono::steady_clock::time_point tpnew = chrono::steady_clock::now();
        chrono::steady_clock::duration duration_in_ites = tpnew - tpold;
        cout << "Time_used_in_iteration " <<  ite << " = " <<  chrono::duration_cast<chrono::seconds>(duration_in_ites).count() << endl;
    } //end the iterations
    
    /* option part ,
     output all the potentials and densities */
//     for(int i = 0 ;i < max_L; i++){
//         ofstream outfile;
//         ofstream outfile2;
//         outfile.open("density" + to_string(i) + ".txt");
//         outfile2.open("potential" + to_string(i) + ".txt");
//         for(int j = 0; j < N; j++){
//             outfile << fx[j] << ' ' << dens[i][j] << ' ' << denv[i][j] << ' ' << den3[i][j] << ' ' << denp[i][j] << endl;
//             outfile2 << fx[j] << ' ' << Phi[i][j] << ' ' << W[i][j] << ' ' << B[i][j] << ' ' << A[i][j] << endl;
//         }
//         outfile.close();
//         outfile2.close();
//     }
//   if (max_L > 1){
//       ofstream outfile;
//       outfile.open("output/density" + to_string(Deformation_parameter) + ".txt");
//       for(int j = 0; j < N; j++){
//           outfile << fx[j] << ' ' << dens[1][j] << ' ' << denv[1][j] << ' ' << den3[1][j] << ' ' << denp[1][j] << endl;
//       }
//   }

 
    
    
    
    
 /* output wavefunctions */

//     for(int ii = 0; ii < occp.size(); ii ++){
//         eig2 test = occp[ii];
//         Solution temp = Solution(test);
//         temp.get_all_wavefunction();
//         cout << ii << " energy for occp = " << temp.energy << endl;
// //        mkdir("wfp" + to_string(ii));
//         string some = "mkdir " + string("testout/wfp") + to_string(ii);
//         const char * path = some.c_str();
//         system(path);
//         for(int i = 0; i < temp.wavefunctions.size(); i ++){
//             ofstream outfile;
//             outfile.open("testout/wfp" + to_string(ii) + '/' + "wavefunction" + to_string(temp.wavefunctions[i].kappa) + ".txt");
//             for(int j = 0; j < N; j++){
//                 outfile << fx[j] << ' ' << temp.wavefunctions[i].upper[j] << ' ' << temp.wavefunctions[i].lower[j] << ' ' << temp.m << ' ' << temp.energy << endl;
//             }
//             outfile.close();
//         }
//     }
    
    
    // for(int ii = 0; ii < occn.size(); ii ++){
    //     eig2 test = occn[ii];
    //     Solution temp = Solution(test);
    //     temp.get_all_wavefunction();
    //     cout << ii << " energy for occn = " << temp.energy << endl;
    //     string some = "mkdir " + string("testout/wfn") + to_string(ii);
    //     const char * path = some.c_str();
    //     system(path);
    //     for(int i = 0; i < temp.wavefunctions.size(); i ++){
    //         ofstream outfile;
    //         outfile.open("testout/wfn" + to_string(ii) + '/' + "wavefunction" + to_string(temp.wavefunctions[i].kappa) + ".txt");
    //         for(int j = 0; j < N; j++){
    //             outfile << fx[j] << ' ' << temp.wavefunctions[i].upper[j] << ' ' << temp.wavefunctions[i].lower[j] << ' ' << temp.m << ' ' << temp.energy << endl;
    //         }
    //         outfile.close();
    //     }
    // }
    
    
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









