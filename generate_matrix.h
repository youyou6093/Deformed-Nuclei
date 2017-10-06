//
//  generate_matrix.h
//  deform_c++
//
//  Created by Junjie Yang on 4/18/17.
//  Copyright © 2017 Junjie Yang. All rights reserved.
//tested

#ifndef generate_matrix_h
#define generate_matrix_h
#include "simps.h"

extern int N;
extern int max_L;
extern vector<double> fx;
extern double hbarc,b,my_tolerance;
extern unordered_map<string, vector<vector<double>> > States_map;
extern unordered_map<string, double> Energy_map;
extern unordered_map<string, double > Angular_map;
extern const double PI;


//not sure whether this is right
string generate_key(double m,int k1,int k2,int L){
    return to_string(int(m))+".5"+"."+to_string(k1)+"."+to_string(k2)+"."+to_string(L);
}


struct params{
    int n1,n2,k1,k2,s1,s2;
    vector<double>g1,g2,f1,f2;
    double m,diag_energy;
    vector<vector<double>> scalar_p,vector_p;
};   // a struct that stores the params




//calculate the matrix element for a specific L channel
double calculate_matrix_element(struct params & my_params,int L){
    double value1,value2,angular_term,value3;
    double E1,E2,E3;
    vector<double> y1,y2,y3,y4,y5;
    //the spherical symmetric potential part
    if (L==0){    //the original part
        if((my_params.n1==my_params.n2)&&(my_params.s1==my_params.s2)&&(my_params.k1==my_params.k2))
            E1=my_params.diag_energy;
        else
            E1=0.0;   //diagnoal part
        //cout<<"test "<<E1<<endl;
        //cout<<"something3"<<endl;
        //cout<<my_params.scalar_p[L].size()<<endl;
        //cout<<my_params.g1.size()<<endl;
        //cout<<my_params.n1<<' '<<my_params.k1<<' '<<my_params.s1<<endl;
        if(my_params.k1==my_params.k2){
            for(int i=0;i<N;i++){
                //cout<<my_params.scalar_p[L][i]<<' '<<my_params.vector_p[L][i]<<endl;
                //cout<<my_params.g1[i]<<' '<<my_params.g2[i]<<endl;
                //cout<<my_params.scalar_p[L][i]<<' '<<my_params.g1[i]<<endl;
                value1=(my_params.scalar_p[L][i]+my_params.vector_p[L][i])*my_params.g1[i]*my_params.g2[i]/hbarc;
                value2=(-my_params.scalar_p[L][i]+my_params.vector_p[L][i])*my_params.f1[i]*my_params.f2[i]/hbarc;
                //cout<<value1<<' '<<value2<<endl;
                value3=-(1/pow(b,2))*(my_params.g1[i]*my_params.f2[i]+my_params.g2[i]*my_params.f1[i])*fx[i];
                //cout<<value1<<' '<<value2<<' '<<value3<<endl;
                y3.push_back(value1+value2+value3);
            }
            //for(int i=0;i<y3.size();i++)
            //cout<<i<<' '<<fx[i]<<' '<<y3[i]<<endl;
            //for(int i=0;i<N;i++) cout<<y3[i]<<endl;
            //E2=my_spline(y3,fx,my_tolerance).integral(); //calculate the integral using gsl;
            E2=simps(y3,fx);
            //cout<<"compare"<<endl;
            //cout<<E2<<' '<<simps(y3, fx)<<endl;
            y3.clear();
        }
        else
            E2=0.0;
        
        //cout<<E1<<' '<<E2<<endl;
        
        return E1+E2;
    }
    
    
    //the deformation part
    else{
        //cout<<my_params.m<<' '<<my_params.k1<<' '<<my_params.k2<<' '<< L<<endl;
        //cout<<generate_key(my_params.m, my_params.k1, my_params.k2, L)<<endl;
        //cout<<generate_key(my_params.m, my_params.k1, my_params.k2, L)<<endl;
        angular_term=Angular_map.find(generate_key(my_params.m, my_params.k1, my_params.k2, L))->second;
        if(angular_term==0) E3=0.0;
        else{
            for(int i=0;i<N;i++){
                value1=(my_params.scalar_p[L][i]+my_params.vector_p[L][i])*my_params.g1[i]*my_params.g2[i]/hbarc;
                value2=(-my_params.scalar_p[L][i]+my_params.vector_p[L][i])*my_params.f1[i]*my_params.f2[i]/hbarc;
                y4.push_back(value1+value2);
            }
            //E3=my_spline(y4,fx,my_tolerance).integral();
            E3=simps(y4,fx);
            y4.clear();
        }
        return E3*sqrt(4*PI/(2*L+1));
    }
}



//this function,get M[i,j] where M is the final matrix;
//I need to find a smarter way to find all the potentials I need;
void generate_matrix(vector<vector<double>> &M,vector<vector<double>> &scalar_potential,vector<vector<double>> & vector_potential,vector<State> &states,double m,int i,int j){
    unordered_map<string, double>::iterator Angular_ptr;
    unordered_map<string, vector<vector<double>>>::iterator States_ptr;
    unordered_map<string, double>::iterator Energy_ptr;
    
    
    
    //vector<double> y1,y2,y3,y4,y5;
    //int n1,k1,n2,k2,s1,s2;
    //double E1,E2,E3;
    //double value1,value2,value3;
    
    //int size_of_L=int(vector_potential.size());     //I need to modify this part
    double matrix_element=0.0;
    //vector<double> upper1,upper2,lower1,lower2;
    //prepare the wavefunction and the state quantum number;
    //finish preparing;
    struct params my_params;
    
    my_params.m=m;
    my_params.n1=states[i].n;
    my_params.n2=states[j].n;
    my_params.k1=states[i].k;
    my_params.k2=states[j].k;
    my_params.s1=states[i].sign;
    my_params.s2=states[j].sign;
    //cout<<states[i].n<<' '<<states[i].k<<' '<<states[i].sign<<endl;
    //cout<<states[j].key<<endl;
    my_params.diag_energy=Energy_map.find(states[i].key)->second; //the energy for states[i]
    
    States_ptr=States_map.find(states[i].key);
    my_params.g1=States_ptr->second[0];
    my_params.f1=States_ptr->second[1];
    States_ptr=States_map.find(states[j].key);
    my_params.g2=States_ptr->second[0];
    my_params.f2=States_ptr->second[1];
    my_params.scalar_p=scalar_potential;
    my_params.vector_p=vector_potential;
    
    
    
    //cout<<"something1"<<endl;
    
    for(int i=0;i<max_L;i++){
        //cout<<i<<' '<<calculate_matrix_element(my_params, i)<<endl;
        matrix_element+=calculate_matrix_element(my_params, i);
        //get the total matrix element for this channel
    }
    
    M[i][j]=matrix_element;
}

//generate the full matrix for a specific m and a specific particle type(this is determined by the scalar or vector potential....
vector<vector<double>> generate_full_matrix(vector<vector<double>> &scalar_potential,vector<vector<double>> & vector_potential,vector<State> &states,double m){
    int matrix_size=int(states.size());
    vector<vector<double>> M(matrix_size,vector<double>(matrix_size,0.0));
    //use the symmetric of this matrix
    for(int i=0;i<states.size();i++)
        for(int j=i;j<states.size();j++)   //this is half of the matrix
            generate_matrix(M, scalar_potential, vector_potential, states, m, i, j);
    for(int j=0;j<states.size();j++)
        for(int i=0;i<j;i++)
            M[j][i]=M[i][j];               //the other half has already been computed
    return M;
}

#endif /* generate_matrix_h */
