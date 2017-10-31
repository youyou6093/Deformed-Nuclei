//just a copy but it has a lot of useful test functions
//try to construct basis
//b=2.4,rmin=0,rmax=20,N=501

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

using namespace std;

//very beginning setup,set up some global variables
double my_tolerance=1e-4;   //I need to find an alternate integrator
double hbarc=197.326;
double b=2.4,rmin=0.0,rmax=20.0;
int N=501;
int max_L=3;
const double PI=atan(1.0)*4;
double h=(rmax-rmin)/(N-1);
vector<double> start ={450.0,400.0,5.0,0.0}; //Phi W B A
vector<double> fx;
unordered_map<string, vector<vector<double>> > States_map;
unordered_map<string, double> Energy_map;
unordered_map<string, double > Angular_map;


// store the upper part and lower part for a wavefunction
void preprocessing(){
    vector<vector<double>> wave_function;
    vector<double> energy;
    for(int i=0;i<N;i++)
        fx.push_back(h*i+rmin);   //construct the x axis for all the basis
    //constructed the States hash table
    vector<State>states_1 = generate_statesm(0.5,max_L);
    for(int i=0;i<states_1.size();i++){
        dirac_oscillator ho(states_1[i].sign,states_1[i].n,states_1[i].k,b,rmin,rmax,N);
        wave_function.clear();   //empty the array
        wave_function.push_back(ho.upper);
        wave_function.push_back(ho.lower);
        States_map.insert(make_pair(states_1[i].key,wave_function));
        Energy_map.insert(make_pair(states_1[i].key, ho.E));
    }
    //now lets constructed the angular part table
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
// maybe I need to put zero into other L channels,no matter what?
void preprocessing_2(vector<vector<double>> & Phi,vector<vector<double>> &W,vector<vector<double>> &B,vector<vector<double>> & A){
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
    
    vector<double> empty;
    for(int i=0;i<N;i++)
        empty.push_back(0.0);
    
    for(int i=0;i<max_L-1;i++){
        Phi.push_back(empty);
        W.push_back(empty);
        B.push_back(empty);
        A.push_back(empty);
    }
}

vector<double> flat_matrix(vector<vector<double>> &M){
    vector<double> flat_matrix;
    int size=int(M.size());
    cout<<size<<endl;
    for(int i=0;i<size;i++)
        for(int j=0;j<size;j++)
            flat_matrix.push_back(M[i][j]);
    cout<<flat_matrix.size();
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
        scalar_p[L][i]+=Potential[i];
    }
}









int main(){
    
    unordered_map<string, double>::iterator Angular_ptr;
    unordered_map<string, vector<vector<double>>>::iterator States_ptr;
    vector<vector<double>> Phi,W,B,A;
    vector<vector<double>> scalar_p,vector_p,scalar_n,vector_n;
    vector<double> Potential,flat;
    vector<State> States_m;
    vector<vector<double>> M;
    matrix_diag diag = matrix_diag();
    vector<eig>  occp,occn; //occupied states of protons and neutrons
    vector<eig>  occp_raw,occn_raw;
    preprocessing();     //create the states_map and angular_map
    preprocessing_2(Phi, W, B, A);      //initialize the PHI W B A potential
    for(int i=0;i<N;i++)              //initialize the pertubation term
        Potential.push_back(fx[i]*0.0);
    int Potential_channel=1;             //set potential to dipole
    cout<<"finishing processing"<<endl;
    
    //------------------------------------------------------body of my program
    //int ite=0;
    
    
    for(int ite=0;ite<1;ite++){
        //break;
        occp.clear();
        occn.clear();            //make sure they are empty
        //generate the potentials for neutron
        generate_potential(Phi, W, B, A, Potential, scalar_n, vector_n, Potential_channel, 0);
        generate_potential(Phi, W, B, A, Potential, scalar_p, vector_p, Potential_channel, 1);
        for(double m=-4.5;m<4.5;m++){
            States_m=generate_statesm(m, max_L);
            //get the matrix for the specific m,neutron
            M=generate_full_matrix(scalar_n, vector_n, States_m, m);
            flat=flat_matrix(M);
            
            
            //another test
            //ofstream matrix_file;
            //matrix_file.open("matrix.txt");
            //for(int i=0;i<States_m.size();i++){
            //    for(int j=0;j<States_m.size();j++){
            //        matrix_file<<M[i][j]<<' ';
            //    }
            //    matrix_file<<endl;
            //}
            //matrix_file.close();
            //break;
            //finish
            diag=matrix_diag(M,int(States_m.size()));    //diagnolize the matrix, add all the eigenvalues and eigenvectors into M_matrix
            diag.get_results();       //diag.results is a matrix full of eigvalues and eigenvectors.
            occn.insert(occn.end(), diag.results.begin(), diag.results.end());
            //get the matrix for the specific m,proton
            M=generate_full_matrix(scalar_p, vector_p, States_m, m);
            diag=matrix_diag(M,int(States_m.size()));
            diag.get_results();
            occp.insert(occp.end(),diag.results.begin(),diag.results.end());
            //diagnolize the matrix, add all the eigenvalues and eigenvectors into M_matrix
            //sort the two solution  based on the eigenvalues and the get biggest N and Z
        }
    }
    
    
    
    
    cout<<"test"<<endl;
    //test my matrix implementation for m=0.5, neutron
    vector<State> states_test=generate_statesm(0.5,max_L);
    //int matrix_size=int(states_test.size());
    //vector<vector<double>> M(matrix_size,vector<double>(matrix_size,0.0));   //construct a matrix that initialized to zero
    //double M_2[matrix_size*matrix_size];    //use 1-d array instead of matrix
    
    //generate_potential(Phi, W, B, A, Potential, scalar_p, vector_p, Potential_channel,0);
    
    //cout<<scalar_p.size()<<endl;
    M=generate_full_matrix(scalar_n, vector_n, states_test, 0.5);
    //cout<<"finish"<<endl;
    //cout<<vector_p.size()<<endl;
    //cout<<states_test.size()<<endl;
    //for(int i=0;i<N;i++) cout<<scalar_p[0][i]<<' '<<vector_p[0][i]<<endl;
    //cout<<states_test[115].n<<' '<<states_test[115].k<<' '<<states_test[115].sign<<endl;
    //generate_matrix(M, scalar_p, vector_p, states_test, 0.5, 10, 10);
    //generate_matrix(M, scalar_p, vector_p, states_test, 0.5, 115, 115);
    //generate_matrix(M, scalar_p, vector_p, states_test, 0.5, 0, 1);
    cout<<M[0][1]<<endl;
    cout<<M[10][10]<<' '<<M[115][115]<<endl;
//--------------------------------------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    //test my result
    vector<string> keys_1;
    for(auto kv : Angular_map) keys_1.push_back(kv.first);
    vector<string> keys_2;
    for(auto kv: States_map) keys_2.push_back(kv.first);
    
    for(int i=0;i<keys_1.size();i++){
        Angular_ptr=Angular_map.find(keys_1[i]);
        //cout<<keys_1[i]<<' '<<Angular_ptr->second<<endl;
    }
    Angular_ptr=Angular_map.find("0.5.1.1.1");
    //cout<<Angular_ptr->second<<endl;
    //dirac_oscillator test=dirac_oscillator(1,1, -3, b, rmin, rmax, 501);
    //States_ptr=States_map.find("11.-3.1");
    //cout<<(States_ptr->second)[0][55]<<' '<<(States_ptr->second)[1][55]<<endl;
    //for(int i=0;i<N;i++) cout<<test.upper[i]<<endl;
    
    //cout<<test.E<<endl;
    //--------------------------------------------------------------------------------------------------------
    double  x[1000], y[1000];
    for (int i = 0; i < 1000; i++){
        x[i]=i*0.01;
        y[i]=x[i]*cos(x[i]);
    }
    my_spline test2=my_spline(y,x,1000,1e-4);
    //cout<<"test spline:"<<endl;
    //cout<<test.eval(0.35)<<' '<<0.35*cos(0.35)<<endl;
    
    
    //cout<<"test integration:"<<endl;
    //cout<<test.integral(0,9.99)<<endl;
    //cout<<test.integral()<<endl;
    //finished integrator part--------------------------------------------------------------------------------
}




//another test
//vector<string> test;
//for(auto kv : States_map) {
//    test.push_back(kv.first);
//}
//cout<<test[2]<<endl;
//wave_function=States_map.find(test[2])->second;
//dirac_oscillator test2(-1,9,-10,b,rmin,rmax,N);
//for(int i=0;i<wave_function[0].size();i+=10){
//    cout<<wave_function[0][i]<<' '<<test2.upper[i]<<endl;
//}








//for test purpose
//dirac_oscillator test(1,2,1,b,rmin,rmax,N);
//for(int i=0;i<N;i++){
//    cout<<fx[i]<<' '<<test.upper[i]<<' '<<test.lower[i]<<endl;
//}
//cout<<test.E<<endl;
