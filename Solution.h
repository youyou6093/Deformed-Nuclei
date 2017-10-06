//
//  Solution.h
//  deform_c++
//
//  Created by Junjie Yang on 6/18/17.
//  Copyright © 2017 Junjie Yang. All rights reserved.
//

//I need to check whether this is right
//most of the wavefunction for L=0 channel is already checked
#ifndef Solution_h
#define Solution_h
#include "symm.h"
#include "generate_state.h"
#include<vector>
#include "states.cpp"

extern int max_L;
extern int N;
extern unordered_map<string, vector<vector<double>> > States_map;
//extern unordered_map<string, double> Energy_map;
//extern unordered_map<string, double > Angular_map;
using std::vector;
using namespace std;

struct eig2{
    eig solution;
    double m;
};




struct coef_pair{             //takes coefficient and state
    double coefs;
    State state;
};

struct solution_wave_function{                // a struct that stores wave function for one kappa
    int kappa;
    vector<double>upper;
    vector<double>lower;
};

class Solution{
public:
    vector<coef_pair> my_pair;
    double m;     //the quantum number of this solution
    double energy;
    vector<int> kappas;
    State primary_state;
    vector<solution_wave_function> wavefunctions;
    
    
    bool add_kappa(int k){                         //decide whether add new kappa into collection
        for(int i=0;i<kappas.size();i++){
            if (k==kappas[i])
                return true;
        }
        kappas.push_back(k);
        return false;
    }
    
    
    Solution(){
        m=0.0;
        energy=0.0;
    }
    
    
    
    Solution(eig2 result){                //1 eigenvalue of diagnolization matrix
        this->m=result.m;                                //m is the quantum numbrer
        energy=result.solution.eigen_values;               //energy for the solution
        coef_pair temp;                           //holder
        vector<State> states=generate_statesm(m,max_L);             //the states for this m
        for(int i=0;i<states.size();i++){
            if (abs(result.solution.eigen_vectors[i])>1e-5){   //automaticlly avoid small coefs
                //dont use this right now
            }
            temp.coefs=result.solution.eigen_vectors[i];
            temp.state=states[i];
            my_pair.push_back(temp);
            add_kappa(states[i].k);                    //add new kappa into the collection
        }            //now I get the pair
        
        
        
    }
    
    
    void get_primary_state(){
        vector<State> states=generate_statesm(m, max_L);
        int max=0;
        for(int i=0;i<my_pair.size();i++)
            if (abs(my_pair[i].coefs)>abs(my_pair[max].coefs))
                max=i;
        
        
        primary_state=states[max];
        
        
    }
    
    
    solution_wave_function get_wave_function(int kappa){
        vector<double> upper;
        vector<double> lower,g,f;
        unordered_map<string, vector<vector<double>>>::iterator States_ptr;
        coef_pair temp;
        for(int i=0;i<N;i++){
            upper.push_back(0.0);
            lower.push_back(0.0);
        }
        for(int i=0;i<my_pair.size();i++){
            temp=my_pair[i];               //get one element from the pair
            if(temp.state.k==kappa){       //if the kappa of that state is what I want
                States_ptr=States_map.find(temp.state.key);
                g=States_ptr->second[0];   //get the wavefunctions
                f=States_ptr->second[1];
                for(int i=0;i<N;i++){
                    upper[i]+=g[i]*temp.coefs;   //add it to our wavefunction
                    lower[i]+=f[i]*temp.coefs;
                }
            }
        }
        solution_wave_function wavefunction1;
        wavefunction1.kappa=kappa;
        wavefunction1.upper=upper;
        wavefunction1.lower=lower;
        return wavefunction1;
    }
    
    void get_all_wavefunction(){               //get all wavefunctions that I need
        for(int i=0;i<kappas.size();i++){
            wavefunctions.push_back(get_wave_function(kappas[i]));
        }
    }
    
    
    
    
    
    
    
};


#endif /* Solution_h */
