
//try to construct basis
//b=2.4,rmin=0,rmax=20,N=501
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

using namespace std;

//very beginning setup,set up some global variables
const double PI=atan(1.0)*4;//
double my_tolerance=1e-4;   //currently the arccuarrcy of the integrator is not very good

//these are model parameters
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
//----------------------------------------------------

double b=2.4,rmin=0.0,rmax=20.0;         
int N=501;
int max_L=1;                                //might need to change
double h=(rmax-rmin)/(N-1);
int proton_number=20;
int neutron_number=28;
vector<double> start ={450.0,400.0,5.0,0.0}; //initial value of Phi W B A in the wood-saxon form
vector<double> fx;
vector<double> empty;
unordered_map<string, vector<vector<double>> > States_map;
unordered_map<string, double> Energy_map;
unordered_map<string, double > Angular_map;


//construct three unordered map: States, energy and Angular
void preprocessing(){
    vector<vector<double>> wave_function;  //store the upper_part and lower part
    vector<double> energy;                 //stor all the energys
    for(int i=0;i<N;i++)
        fx.push_back(h*i+rmin);   //construct the x array for all the basis
    
    //constructed the States hash table
    vector<State>states_1 = generate_statesm(0.5,max_L);  //this is the largest basis
    for(int i=0;i<states_1.size();i++){
        dirac_oscillator ho(states_1[i].sign,states_1[i].n,states_1[i].k,b,rmin,rmax,N);
        wave_function.clear();   //empty the array
        wave_function.push_back(ho.upper);
        wave_function.push_back(ho.lower);
        States_map.insert(make_pair(states_1[i].key,wave_function));
        Energy_map.insert(make_pair(states_1[i].key, ho.E));
    }
    //now lets constructed the angular part table(the values are read from file)
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


//initialize the potentials
void preprocessing_2(vector<vector<double>> & Phi,vector<vector<double>> &W,vector<vector<double>> &B,vector<vector<double>> & A){
    //double PI=3.0;
    double value;
    Phi.clear();
    W.clear();
    B.clear();
    A.clear();
    //make sure that the potentials are empty at the start of each iteration
    vector<double> phi,w,b,a;
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
    //just try to be safe
    //vector<double> empty;
    for(int i=0;i<N;i++)
        empty.push_back(0.0);  //now empty is a global variable
    
    for(int i=0;i<max_L-1;i++){
        Phi.push_back(empty);
        W.push_back(empty);
        B.push_back(empty);
        A.push_back(empty);
    }
}


//flat a matrix into an vector
vector<double> flat_matrix(vector<vector<double>> &M){
    vector<double> flat_matrix;
    int size=int(M.size());
    
    for(int i=0;i<size;i++)
        for(int j=0;j<size;j++)
            flat_matrix.push_back(M[i][j]);
    if(flat_matrix.size()!=size*size) cout<<"error"<<endl;   //try to be safe
    return flat_matrix;
}

    




//generate the scalar and vector potentials
void generate_potential(vector<vector<double>> &Phi,vector<vector<double>> &W, vector<vector<double>> &B,vector<vector<double>> &A,vector<double> &Potential,vector<vector<double>> &scalar_p,vector<vector<double>> &vector_p,int L,int particle_type){//L represent the potential channel
    
    scalar_p.clear();
    vector_p.clear();
    //clear the potential so that it will be empty at the start of each iterations
    
    vector<double> temp1,temp2;
    if (particle_type==0){                                //neutron
        for(int i=0;i<Phi.size();i++){                    //for all the channel
            //temp1 and temp2 is the potential for 1 channel
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
    
    for(int i=0;i<N;i++){
        scalar_p[L][i]+=Potential[i];                    //put the excitation potential into the scarlar potential
    }
}




bool compare_eig2(eig2 a,eig2 b){
    return (a.solution.eigen_values<b.solution.eigen_values);
}

void get_solution(vector<eig2> &occp_raw,vector<eig2> &occn_raw,vector<eig2> &occp,vector<eig2> &occn){
    //filter the solution,only accept that are physical
    for(int i=0;i<occp_raw.size();i++){
        occp_raw[i].solution.eigen_values=occp_raw[i].solution.eigen_values*hbarc-939;
        if ((occp_raw[i].solution.eigen_values<0) && (occp_raw[i].solution.eigen_values>-939)){
            occp.push_back(occp_raw[i]);
        }
        
    }
    for(int i=0;i<occn_raw.size();i++){
        occn_raw[i].solution.eigen_values=occn_raw[i].solution.eigen_values*hbarc-939;
        if ((occn_raw[i].solution.eigen_values<0) && (occn_raw[i].solution.eigen_values>-939)){
            occn.push_back(occn_raw[i]);
        }
    }
    
    //sort the solution, and get the states that are occupied
    sort(occn.begin(),occn.end(),compare_eig2);
    sort(occp.begin(),occp.end(),compare_eig2);
    if(occn.size()>neutron_number) occn.erase(occn.begin()+neutron_number,occn.end());
    if(occp.size()>proton_number) occp.erase(occp.begin()+proton_number,occp.end());
    
}




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


//vector<Solution> get_solutions_object(vector<eig2> occ){
//    vector<Solution> Final_occ;
//    Solution temp;
//    for(int i=0;i<occ.size();i++){
//        temp=Solution(occ[i]);
//        Final_occ.push_back(temp);
//    }
//    return Final_occ;
//}



int main(){
    //test the time
    chrono::steady_clock::time_point tp1 = chrono::steady_clock::now();
    unordered_map<string, double>::iterator Angular_ptr;
    unordered_map<string, vector<vector<double>>>::iterator States_ptr;
    vector<vector<double>> Phi,W,B,A,dens,denv,denp,den3,EFF_Phi,EFF_W,EFF_B,EFF_A;
    vector<vector<double>> scalar_p,vector_p,scalar_n,vector_n;
    vector<double> Potential,flat; //flat is the flatted matrix
    vector<State> States_m;        //the basis quantum number of a specific m
    vector<vector<double>> M;      //Matrix to be diagnolized
    matrix_diag diag = matrix_diag();      //diagnolization class
    vector<eig2>  occp,occn; //occupied states of protons and neutrons
    vector<eig2>  occp_raw,occn_raw;      //raw solution of matrix
    vector<eig2> temp_solution;
    //determine the Z and N;
    cout << "Z=: ";
    cin  >> proton_number;
    cout << "N=: ";
    cin >> neutron_number;
    //finish
    preprocessing();                    //create the states_map and angular_map
    preprocessing_2(Phi, W, B, A);      //initialize the PHI W B A potential
    for(int i=0;i<N;i++)                //initialize the excitation term
        Potential.push_back(fx[i]*0.001);                                        //might need to change
    int Potential_channel=0;            //set potential to dipole,might need to change
    //get empty density;
    for(int i=0;i<max_L;i++){
        dens.push_back(empty);
        denv.push_back(empty);
        denp.push_back(empty);
        den3.push_back(empty);
    }
    EFF_Phi=dens; EFF_W=dens;EFF_A=dens;EFF_B=dens;   //set everything to zero;
    //determine the range of m;
    double min_m = -magic(max(proton_number,neutron_number));
    double max_m = -min_m;
    cout<<min_m<<endl;
    
    //------------------------------------
    cout<<"finishing processing"<<endl;
    
    //------------------------------------------------------body of my program
    
    vector<Solution> Final_occp,Final_occn;
    
    for(int ite=0;ite<50;ite++){
        occp_raw.clear();
        occn_raw.clear();
        occp.clear();
        occn.clear();            //make sure they are empty
        //generate the potentials for neutron
        generate_potential(Phi, W, B, A, Potential, scalar_n, vector_n, Potential_channel, 0);
        generate_potential(Phi, W, B, A, Potential, scalar_p, vector_p, Potential_channel, 1);
        for(double m = min_m ; m < max_m + 1 ; m++){
            //cout<<"m="<<m<<endl;
            States_m=generate_statesm(m, max_L);
            M=generate_full_matrix(scalar_n, vector_n, States_m, m);  //get the matrix for the specific m,neutron
            if (M.size()!=States_m.size()) cout<<"error!!"<<endl;   //try to be safe
            flat=flat_matrix(M);
            diag=matrix_diag(flat,int(States_m.size()));    //diagnolize the matrix, add all the eigenvalues and eigenvectors into M_matrix
            diag.get_results();       //diag.results is a matrix full of eigvalues and eigenvectors.
            temp_solution=get_temp_solution(diag.results, m);
            occn_raw.insert(occn_raw.end(), temp_solution.begin(), temp_solution.end());
            M=generate_full_matrix(scalar_p, vector_p, States_m, m);    //get the matrix for the specific m,proton
            flat=flat_matrix(M);
            diag=matrix_diag(flat,int(States_m.size()));
            diag.get_results();
            temp_solution=get_temp_solution(diag.results, m);
            occp_raw.insert(occp_raw.end(),temp_solution.begin(),temp_solution.end());
            //diagnolize the matrix, add all the eigenvalues and eigenvectors into M_matrix
            //sort the two solution  based on the eigenvalues and the get biggest N and Z
        }
        get_solution(occp_raw, occn_raw, occp, occn);
        //for(int i=0;i<occp.size();i++) cout<<occp[i].solution.eigen_values<<' '<<occp[i].m<<endl;
        //for(int i=0;i<occn.size();i++) cout<<occn[i].solution.eigen_values<<' '<<occn[i].m<<endl;
        Final_occn=get_solutions_object(occn);
        Final_occp=get_solutions_object(occp);
        //next step: test whether it reproduce the right wave functions
        //for(int i=0;i<occp.size();i++) cout<<Final_occp[i].energy<<' '<<Final_occp[i].m<<endl;
        //for(int i=0;i<occn.size();i++) cout<<Final_occn[i].energy<<' '<<Final_occn[i].m<<endl;
        //structure of my next part
        //vector<Solution> solution = get_solutions_object(occp);
        
//        Final_occp[0].get_all_wavefunction();
//        for(int i=0;i<Final_occp[0].kappas.size();i++)  cout<<Final_occp[0].kappas[i]<<endl;
//        Final_occp[0].get_primary_state();
//        cout<<Final_occp[0].primary_state<<endl;
//        solution_wave_function test1=Final_occp[0].get_wave_function(-1);
//        cout<<test1.upper[100]<<' '<<test1.lower[100]<<endl;
//        for( int i=0;i<N;i++){
//            cout<<test1.upper[i]<<' '<<test1.lower[i]<<endl;
//        }
        
//        Density test=Density(occn);
//        test.Compute_all();
////        for(int i=0;i<N;i++)
////            if(i%50==0)
////                cout<<test.density[0].s[i]<<' '<<test.density[0].v[i]<<endl;
        generate_density(occn,occp,dens,denv,denp,den3); //calculate all the density based on the occupied states
        update_potential(EFF_Phi, EFF_B, EFF_A, EFF_W, Phi, W, B, A, dens, denv, den3, denp);
        cout<<"ite="<<ite<<endl;
        for(int i=0;i<max_L;i++)
            cout <<"channel="<<i<<':'<< Phi[i][0] <<' '<< W[0][0]<<' '<<B[0][0]<<' '<<A[0][0]<<endl;
        cout<<occp.size()<<' '<<occn.size()<<endl;
        cout<<"E/A="<<compute_energy(occp, occn, Phi, W, B, A, dens, denv, den3, denp)<<endl;
        //ofstream myfile;
        //myfile.open("potential.txt");
        //for(int i=0;i<N;i++)
        //    myfile<<fx[i]<<' '<<Phi[0][i]<<' '<<W[0][i]<<' '<<B[0][i]<<' '<<A[0][i]<<endl;
        
        //double test_y[3]={0,1};
        //double test_x[3]={0,1};
        //cout<<my_spline(test_y,test_x,2,my_tolerance).integral()<<endl;
//        ofstream myfile;
//        myfile.open("density.txt");
//        for(int i=0;i<N;i++){
//            myfile<<fx[i]<<' '<<dens[0][i]<<' '<<denv[0][i]<<' '<<den3[0][i]<<' '<<denp[0][i]<<endl;
//        }
//        myfile.close();
        //update_potential(EFF_Phi,EFF_B,EFF_A,EFF_W,Phi,W,B,A,dens,denv,den3,denp);
        
//        test for the three j symbol
//        double a1=1.0,a2=2.0,a3=3.0;
//        double b1=0.,b2=0.,b3=0.;
//        
//        double result = -1.0;
//        result = tj_(&a1,&a2,&a3,&b1,&b2,&b3);//three J symbol
//        std::cout << result << std::endl;
//        
//        
//        double j1=0.5,j2=0.5,j3=0.5,j12=1.,j23=1.0,j=1.5;
//        std::cout << sj_(&j1,&j2,&j12,&j3,&j,&j23) << std::endl;//six j symbol
//        double a4 = 0.5,a5=0.5,a6=0.0,a7=1.5,a8=2.5,a9=3.0;
//        std::cout << coef9_(&a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8,&a9) << std::endl;//nine j symbol
//        return 0;

    }
    chrono::steady_clock::time_point tp2 = chrono::steady_clock::now();
    chrono::steady_clock::duration d = tp2-tp1;
    cout <<"Total time used:"<<chrono::duration_cast<chrono::seconds>(d).count()<<endl;
        
        
        
//        //test part
//        Final_occn[6].get_primary_state();
//        //cout<<Final_occn[0].primary_state<<endl;
//        solution_wave_function test=Final_occn[6].get_wave_function(Final_occn[6].primary_state.k);
//        ofstream myfile;
//        myfile.open("output.txt");
//        for(int i=0;i<test.upper.size();i++){
//            //cout<<fx[i]<<' '<<test.upper[i]<<' '<<test.lower[i]<<endl;
//            myfile<<fx[i]<<' '<<test.upper[i]<<' '<<test.lower[i]<<endl;
//        }
//        myfile.close();
        
    
    
    
    
    //cout<<endl;
    //cout<<"test"<<endl;
    //test my matrix implementation for m=0.5, neutron
    //vector<State> states_test=generate_statesm(0.5,max_L);
    //M=generate_full_matrix(scalar_n, vector_n, states_test, 0.5);
    //cout<<M[0][1]<<endl;
    //cout<<M[10][10]<<' '<<M[115][115]<<endl;
    //--------------------------------------------------------------------------------------------------------
}


