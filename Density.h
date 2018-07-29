//
//  Density.h
//  deform_c++
//
//  Created by Junjie Yang on 7/13/17.
//  Copyright © 2017 Junjie Yang. All rights reserved.
#ifndef Density_h
#define Density_h
#include <math.h>
#include "Solution.h"
#include <vector>
#include "vector_compute.h"
#include "generate_matrix.h"



/*This function transfer the eigenvalues and eigenvetors into an object named solution*/
vector<Solution> get_solutions_object(vector<eig2> occ){
    vector<Solution> Final_occ;
    Solution temp;
    for(int i=0;i<occ.size();i++){
        temp=Solution(occ[i]);
        Final_occ.push_back(temp);
    }
    return Final_occ;
}


/*get possilbe L channel for k1 and k2 state
  the parity cg has already been considered
 */
vector<int> get_possible_L(int k1,int k2){
    int j1 = 2 * abs(k1) - 1;
    int j2 = 2 * abs(k2) - 1;
    vector<int> L_group;
    
    /*get l1 and l2*/
    int l1, l2;
    if (k1 > 0) l1 = k1;
    else l1 = -k1 - 1;
    if (k2 > 0) l2 = k2;
    else l2 = -k2 - 1;
    
    /*parity cg coefs*/
    for(int i = abs(j1-j2) ; i <= j1 + j2 ; i+=2){
        if(check(l1, l2, i / 2) == 1)
            L_group.push_back(i / 2);
    }
    return L_group;
}


/* A Struct that stores the vector density and scalar density*/
struct den{
    vector<double> s;
    vector<double> v;
};




class Density{
public:
    vector<Solution> solution;
    vector<den> density;
    
    
    Density(vector<eig2> & occ){
        solution=get_solutions_object(occ);  //Get a set of obeject of class of solution
        for(int i=0;i<solution.size();i++)
            solution[i].get_all_wavefunction(); //Compute all the wavefunction
        vector<double> empty;
        den big_empty;
        for(int i=0;i<N;i++)
            empty.push_back(0.0);
        big_empty.s=empty;
        big_empty.v=empty;
        for(int i=0;i<max_L;i++){
            density.push_back(big_empty);    //Give zero to all density values
        }
    }
    
    
    /*add the new desnity calculated form 1 state*/
    void append(den temp,int channel){
        for(int i=0;i<N;i++){
            density[channel].s[i]+=temp.s[i];
            density[channel].v[i]+=temp.v[i];
        }
    }
    
    /*pick 1 state and 2 kappa*/
    void compute_one_(int a,int b,int num){
        Solution x=solution[num];      //one of the occ state
        int m=x.m;                     //get m
        double coef;             //coef is the A(m,k1,k2,L) in my calculations
        vector<int> kappas=x.kappas;   //all Kappa in this occ state
        int k1=x.wavefunctions[a].kappa;        // the first kappa of this state
        int k2=x.wavefunctions[b].kappa;        // the second kappa of this state
        vector<double> radial1=vector_multiple(x.wavefunctions[a].upper,x.wavefunctions[b].upper);          //f1*f1prime
        vector<double> radial2=vector_multiple(x.wavefunctions[a].lower,x.wavefunctions[b].lower);          //f2*f2prime
        for(int i=0;i<N;i++){
            radial1[i]/= (fx[i] * fx[i]);    //divide by r**2
            radial2[i]/= (fx[i] * fx[i]);
        }
        radial1[0]=radial1[1];
        radial2[0]=radial2[1];
        vector<int> L_group = get_possible_L(k1,k2);  //all the possible L get from this two kappa
        for (int i=0;i<L_group.size();i++){
            /*L_group[i] is the channel that I need to add the denisty*/
            if(L_group[i]<max_L){
                coef = Angular_depedencek(k1, k2, m, L_group[i]);
                den temp;
//                temp.s.clear();
//                temp.v.clear();
                temp.s.reserve(N);
                temp.v.reserve(N);
                for(int j=0;j<N;j++){
                    temp.s.push_back(coef*(radial1[j]-radial2[j]));          //scalar
                    temp.v.push_back(coef*(radial1[j]+radial2[j]));          //vector
                }
                append(temp, L_group[i]);
            }
        }
    }
    
    /*pick one states*/
    void compute(int num){
        Solution x=solution[num];
        for(int i=0;i<x.kappas.size();i++)
            for(int j=0;j<x.kappas.size();j++)
                compute_one_(i, j, num);
    }
    
    

    /*pick all the states*/
    void Compute_all(){
        for(int i=0;i<solution.size();i++){
            compute(i);
        }
    }
};



/*Get density of everything for every iteration*/
void generate_density(vector<eig2> &occn, vector<eig2> &occp, vector<vector<double>> &dens,
                          vector<vector<double>> &denv,vector<vector<double>> &denp,vector<vector<double>> &den3){
    Density all_states_p(occp),all_states_n(occn); // the constructor will find the wave funcitons of all occ state
    vector<den> density_p,density_n;
    vector<vector<double>> denvp,denvn,densn,densp;
    
    all_states_p.Compute_all();
    all_states_n.Compute_all();
    density_p=all_states_p.density;
    density_n=all_states_n.density;

    for(int i=0;i<max_L;i++){
        for(int j=0;j<N;j++){
            dens[i][j]=density_n[i].s[j] + density_p[i].s[j];
            denv[i][j]=density_n[i].v[j] + density_p[i].v[j];
            den3[i][j]=density_p[i].v[j] - density_n[i].v[j];
            denp[i][j]=density_p[i].v[j];
        }
    }

}





#endif /* Density_h */
