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

/*set up global varialbes*/
const double PI=atan(1.0)*4;
string parameters = "parameter.txt";  // the file contains the model parameters
double my_tolerance=1e-4;   //tolerance of the integrator

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
/*Max value of L CHANNEL, it is initialized to 1( nothing at all! ). However, it will change after 50(or other values) iterations */
int max_L = 1;
/*The potential I want to add, initialized to 0 but will change to 1(or other value) after some iterations */
int Potential_channel = 0;
/*After a few iterations it will also change */
double Deformation_parameter = 0.00;
double h=(rmax-rmin)/(N-1);       //grid size
int proton_number, neutron_number;
vector<double> start ={450.0,400.0,5.0,0.0}; //initial value of Phi W B A in the wood-saxon form
vector<double> fx;                //X value array
vector<double> empty;             //empty vector, will be initialized later
unordered_map<string, vector<vector<double>> > States_map;   //store the infos about the wavefunctions when given states
unordered_map<string, double> Energy_map;                    //store the energies for every state
unordered_map<string, double > Angular_map;                  //Stroe the value A in my notes


/*construct three unordered map: States, energy and Angular and get global data */
void preprocessing(){
/* get the input parameters */
    ifstream parameter_file;
    parameter_file.open(parameters);
    parameter_file >> ms;         
    parameter_file >> mv;         
    parameter_file >> mp;         
    parameter_file >> mg;         
    parameter_file >> gs;
    parameter_file >> gv;
    parameter_file >> gp;
    parameter_file >> gg;
    parameter_file >> lambdas;
    parameter_file >> lambdav;
    parameter_file >> lambda;
    parameter_file >> ks;
    parameter_file >> ka;   
/*--------------------------------------------- */

    cout << ms << ' ' << mv << ' ' << mp << ' ' << mg << endl;
    cout << gs << ' ' << gv << ' ' << gp << ' ' << gg << endl;
    cout<< lambdas << ' ' << lambdav << ' ' << lambda << ' ' << ks << ' ' << ka << endl; 
    ms /= hbarc;
    mv /= hbarc;
    mp /= hbarc;
    mg /= hbarc;
    ka /= hbarc;
    
    
    vector<vector<double>> wave_function;  //store the upper_part and lower part
    vector<double> energy;                 //stor all the energys
    for(int i=0;i<N;i++)
        fx.push_back(h*i+rmin);   //construct the x array for all the basis
    
    /*constructed the States hash table*/
    vector<State>states_1 = generate_statesm(0.5,MAX_KAPPA);  //this is the largest basis, gives me everything in the basis
    for(int i=0;i<states_1.size();i++){
        dirac_oscillator ho(states_1[i].sign,states_1[i].n,states_1[i].k,b,rmin,rmax,N);
        wave_function.clear();   //empty the array
        wave_function.push_back(ho.upper);
        wave_function.push_back(ho.lower);
        States_map.insert(make_pair(states_1[i].key, wave_function));
        Energy_map.insert(make_pair(states_1[i].key, ho.E));
    }
    /*construct the angular hash table, basically read values from the file and store it inside unordered set*/
    vector<string> keys_1,value_1;
    string word;
    ifstream infile;
    infile.open("test.txt");
    while(infile>>word){
        keys_1.push_back(word);
        infile>>word;
        value_1.push_back(word);
    }
    infile.close();
    for(int i=0;i<keys_1.size();i++){
        Angular_map.insert(make_pair(keys_1[i],stod(value_1[i])));
    }
}


/*initialize the potentials to wood-saxon form*/
void preprocessing_2(vector<vector<double>> & Phi,vector<vector<double>> &W,vector<vector<double>> &B,vector<vector<double>> & A){
    double value;
    
    /* clear the Potentials */
    Phi.clear();
    W.clear();
    B.clear();
    A.clear();
    vector<double> phi,w,b,a;
    /* wood saxon part */
    for(int i=0;i<N;i++){
        value=(exp(fx[i]-4.3)/0.7+1);
        phi.push_back(start[0]/value);
        w.push_back(start[1]/value);
        b.push_back(start[2]/value);
        a.push_back(start[3]/value);
    }
    Phi.push_back(phi);
    W.push_back(w);
    B.push_back(b);
    A.push_back(a);
    phi.clear();
    w.clear();
    b.clear();
    a.clear();
    /*initialize the 0 channel to wood saxon form*/
    
    /*initialize all other channels to 0*/
    for(int i=0;i<N;i++)
        empty.push_back(0.0);
    for(int i=0;i<max_L-1;i++){         //max_L - 1 because I already assign the 0 channel
        Phi.push_back(empty);
        W.push_back(empty);
        B.push_back(empty);
        A.push_back(empty);
    }
}


/*flat a matrix into an vector*/
vector<double> flat_matrix(vector<vector<double>> &M){
    vector<double> flat_matrix;
    int size=int(M.size());
    for(int i=0;i<size;i++)
        for(int j=0;j<size;j++)
            flat_matrix.push_back(M[i][j]);
    if(flat_matrix.size()!=size*size) cout<<"error"<<endl;   //try to be safe
    return flat_matrix;
}



/*generate the scalar and vector potentials based on PHI W B A
 L is the channel of the potential we add*/
void generate_potential(vector<vector<double>> &Phi,vector<vector<double>> &W, vector<vector<double>> &B,vector<vector<double>> &A,vector<double> &Potential,vector<vector<double>> &scalar_p,vector<vector<double>> &vector_p,int L,int particle_type){
    
    /*clear the scalar potential and vector potential at the beginning, very important*/
    scalar_p.clear();
    vector_p.clear();
    vector<double> temp1,temp2;
    if (particle_type==0){                                //neutron
        for(int i=0;i<Phi.size();i++){                    //for all the channel
            /*temp1 and temp2 is the potential for 1 channel*/
            temp1.clear();                                //clear the container
            temp2.clear();
            for(int j=0;j<Phi[0].size();j++){
                temp1.push_back(-Phi[i][j]);
                temp2.push_back(W[i][j]-0.5*B[i][j]);
            }
            scalar_p.push_back(temp1);
            vector_p.push_back(temp2);
        }
    }
    else{                                                //proton
        for(int i=0;i<Phi.size();i++){
            temp1.clear();
            temp2.clear();
            for(int j=0;j<Phi[0].size();j++){
                temp1.push_back(-Phi[i][j]);
                temp2.push_back(W[i][j]+0.5*B[i][j]+A[i][j]);
            }
            scalar_p.push_back(temp1);
            vector_p.push_back(temp2);
        }
    }
    
    
    /* put the potential we add to the scalar and vector part, depend on what I want the potential to be*/
    for(int i=0;i<N;i++){
        scalar_p[L][i]+=Potential[i];
    }
}



/*compare function for eig2 struct*/
bool compare_eig2(eig2 a,eig2 b){
    return (a.solution.eigen_values<b.solution.eigen_values);
}


/* after we get all the raw solutions, change the solution to MEV and only keep the energies between -939 and 5,
 make sure all the energy is negative in the end*/
void get_solution(vector<eig2> &occp_raw,vector<eig2> &occn_raw,vector<eig2> &occp,vector<eig2> &occn){
    for(int i=0;i<occp_raw.size();i++){
        occp_raw[i].solution.eigen_values=occp_raw[i].solution.eigen_values*hbarc-939;
        if ((occp_raw[i].solution.eigen_values<5) && (occp_raw[i].solution.eigen_values>-939)){
            occp.push_back(occp_raw[i]);
        }
        
    }
    for(int i=0;i<occn_raw.size();i++){
        occn_raw[i].solution.eigen_values=occn_raw[i].solution.eigen_values*hbarc-939;
        if ((occn_raw[i].solution.eigen_values<5) && (occn_raw[i].solution.eigen_values>-939)){
            occn.push_back(occn_raw[i]);
        }
    }
    
    /*sort the solutions*/
    sort(occn.begin(),occn.end(),compare_eig2);
    sort(occp.begin(),occp.end(),compare_eig2);
    /*delete some solutions if I already got enough occ states*/
    if(occn.size()>neutron_number) occn.erase(occn.begin()+neutron_number,occn.end());
    if(occp.size()>proton_number) occp.erase(occp.begin()+proton_number,occp.end());
    
}



/* for specific m, return the solution eig2,contains m */
vector<eig2> get_temp_solution(vector<eig> &results,double m){
    eig2 temp;
    vector<eig2> temp_solution;
    for(int i=0;i<results.size();i++){
        temp.m=m;
        temp.solution=results[i];
        temp_solution.push_back(temp);
    }
    return temp_solution;
}




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
    cout<<min_m<<endl;
    //------------------------------------
    cout<<"finishing processing"<<endl;
   
    
    
    /*------------------------------------------------------body of my program*/
    vector<Solution> Final_occp,Final_occn;                 //Final solution for 1 iteration
    
    
    
    for(int ite=0;ite<2;ite++){
        
        chrono::steady_clock::time_point tpold = chrono::steady_clock::now();
        
        occp_raw.clear();
        occn_raw.clear();
        occp.clear();
        occn.clear();            //make sure they are empty
        /* the newly implemented part, when
         The iteration is bigger thant some presetted time,
         I will add dipole potential and add more L channel*/
        if (ite > -1){
            max_L = 3;
            Deformation_parameter =  0.1 ;
            Potential_channel = 1;
            for (int i = 1 ; i < max_L; i++){
                Phi.push_back(empty);
                W.push_back(empty);
                B.push_back(empty);
                A.push_back(empty);
                dens.push_back(empty);
                denv.push_back(empty);
                denp.push_back(empty);
                den3.push_back(empty);
            }
            EFF_Phi=dens; EFF_W=dens;EFF_A=dens;EFF_B=dens;
        }
        for(int i=0;i<N;i++)                //initialize the excitation term
            Potential[i] = fx[i] * Deformation_parameter;
            //might need to change
        /*----------------------------------*/
        generate_potential(Phi, W, B, A, Potential, scalar_n, vector_n, Potential_channel, 0);   //neutron
        generate_potential(Phi, W, B, A, Potential, scalar_p, vector_p, Potential_channel, 1);   //proton
        
        /* solve the dirac equation for one iteration */
        for(double m = min_m ; m < max_m + 1 ; m++){
            States_m=generate_statesm(m, max_L);
            generate_full_matrix(scalar_n, vector_n, States_m, m);
            vector<vector<double>> M = generate_full_matrix(scalar_n, vector_n, States_m, m);  //get the matrix for the specific m,neutron
            if (M.size()!=States_m.size()) cout<<"error!!"<<endl;   //try to be safe
            /*another test */
            flat=flat_matrix(M);
            diag=matrix_diag(flat,int(States_m.size()));    //diagnolize the matrix, add all the eigenvalues and eigenvectors into M_matrix
            diag.get_results();                             //diag.results is a matrix full of eigvalues and eigenvectors.
            temp_solution=get_temp_solution(diag.results, m);   //solution for this m
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
            outfile << dens[i][j] << ' ' << denv[i][j] << ' ' << den3[i][j] << ' ' << denp[i][j] << endl;
            outfile2 << Phi[i][j] << ' ' << W[i][j] << ' ' << B[i][j] << ' ' << A[i][j] << endl;
        }
        outfile.close();
        outfile2.close();
    }
    
    /* get time */
    chrono::steady_clock::time_point tp2 = chrono::steady_clock::now();
    chrono::steady_clock::duration d = tp2-tp1;
    cout <<"Total time used:"<<chrono::duration_cast<chrono::seconds>(d).count()<<endl;
}


















/* output potential */
//ofstream myfile;
//myfile.open("potential.txt");
//for(int i=0;i<N;i++)
//    myfile<<fx[i]<<' '<<Phi[0][i]<<' '<<W[0][i]<<' '<<B[0][i]<<' '<<A[0][i]<<endl;


/* output density */
//        ofstream myfile;
//        myfile.open("density.txt");
//        for(int i=0;i<N;i++){
//            myfile<<fx[i]<<' '<<dens[0][i]<<' '<<denv[0][i]<<' '<<den3[0][i]<<' '<<denp[0][i]<<endl;
//        }
//        myfile.close();

/* output the single particle energy */
//cout << "proton_energy" << endl;
//for(int i=0;i<occp.size();i++) cout<<Final_occp[i].energy<<' '<<Final_occp[i].m<<endl;
//cout << "neutron_energy" << endl;
//for(int i=0;i<occn.size();i++) cout<<Final_occn[i].energy<<' '<<Final_occn[i].m<<endl;
/* ----------------------------------- */


//        /*test*/
//        Solution test = Final_occp[0];
//        test.get_primary_state();
//        cout << test.energy << ' ' << test.m << ' ' << test.primary_state << endl;
//        test.get_all_wavefunction();
//        for( int i =0; i < test.my_pair.size(); i++)
//            cout << test.my_pair[i].coefs << ' ' << test.my_pair[i].state << endl;
//        for( int i = 0; i < test.wavefunctions.size(); i++){
//            ofstream infile;
//            infile.open(to_string(test.wavefunctions[i].kappa) + "wavefunctionc.txt");
//            for( int x = 0 ; x < N; x++){
//                infile << fx[x] << ' ' << test.wavefunctions[i].upper[x] << ' ' << test.wavefunctions[i].lower[x] << endl;
//            }
//            infile.close();
//        }
////        for(int i = 0; i < N; i++){
////            cout << scalar_n[0][i] << ' ' << vector_n[0][i] << ' ' << scalar_p[0][i] << ' ' << vector_p[0][i] << endl;
////        }
//        break;
//
