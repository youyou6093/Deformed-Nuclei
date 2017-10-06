//
//  Density.h
//  deform_c++
//
//  Created by Junjie Yang on 7/13/17.
//  Copyright © 2017 Junjie Yang. All rights reserved.
// need to check for correctnese
//tested for the naive case

#ifndef Density_h
#define Density_h
#include "Solution.h"
#include <vector>
#include "vector_compute.h"
#include "generate_matrix.h"

vector<Solution> get_solutions_object(vector<eig2> occ){
    vector<Solution> Final_occ;
    Solution temp;
    for(int i=0;i<occ.size();i++){
        temp=Solution(occ[i]);
        Final_occ.push_back(temp);
    }
    return Final_occ;
}

vector<int> get_possible_L(int k1,int k2){
    double j1=abs(k1)-0.5;
    double j2=abs(k2)-0.5;
    vector<int> L_group;
    for(int i=int(abs(j1-j2)) ; i<int(abs(j1+j2+1)) ; i++){
        L_group.push_back(i);
    }
    return L_group;
}

struct den{
    vector<double> s;
    vector<double> v;
    //int channel;
};




class Density{
public:
    vector<Solution> solution;
    vector<den> density;
    
    
    Density(vector<eig2> & occ){
        solution=get_solutions_object(occ);   //get a set of obeject of class of solution
        for(int i=0;i<solution.size();i++)
            solution[i].get_all_wavefunction();           //compute all the wavefunction
        vector<double> empty;
        den big_empty;
        for(int i=0;i<N;i++)
            empty.push_back(0.0);
        big_empty.s=empty;
        big_empty.v=empty;
        for(int i=0;i<max_L;i++){
            density.push_back(big_empty);    //give zero to all density values
        }
    }
    
    
    //add the new desnity calculated form 1 single occupied state
    void append(den temp,int channel){
        for(int i=0;i<N;i++){
            density[channel].s[i]+=temp.s[i];
            density[channel].v[i]+=temp.v[i];
        }
    }
    
    
    void compute_one_(int a,int b,int num){
        Solution x=solution[num];      //one of the occ state
        double m=x.m;                  //get m
        double coef,coef2;             //coef is the A(m,k1,k2,L) in my calculations
        vector<int> kappas=x.kappas;   //all Kappa in this occ state
        int k1=x.wavefunctions[a].kappa;        // the first kappa of this state
        int k2=x.wavefunctions[b].kappa;        // the second kappa of this state
        //cout<<"k1="<<k1<<' '<<"k2="<<k2<<endl;
        //cout<<"pos1"<<endl;
        //cout<<x.wavefunctions[a].upper[100]<<' '<<x.wavefunctions[a].lower[100]<<endl;
        vector<double> radial1=vector_multiple(x.wavefunctions[a].upper,x.wavefunctions[b].upper);          //f1*f1prime
        vector<double> radial2=vector_multiple(x.wavefunctions[a].lower,x.wavefunctions[b].lower);          //f2*f2prime
        for(int i=0;i<N;i++){
            radial1[i]/=pow(fx[i],2);    //divide by r**2
            radial2[i]/=pow(fx[i],2);
        }
        radial1[0]=radial1[1];
        radial2[0]=radial2[1];
        //cout<<"pos2"<<endl;
        //cout<<radial1[100]<<' '<<radial2[100]<<endl;
        vector<int> L_group = get_possible_L(k1,k2);  //all the possible L get from this two kappa
        for (int i=0;i<L_group.size();i++){
            if(L_group[i]<max_L){                             //make sure it is within range
                //cout<<m<<' '<<k1<<' '<<k2<<' '<<L_group[i]<<endl;
                coef = Angular_map.find(generate_key(m, k1, k2, L_group[i]))->second;
                
                coef2 = sqrt((2.0*i+1)/4/PI);
                //cout<<coef2<<endl;
                den temp;                       //temp holder for this channel of this state
                for(int j=0;j<N;j++){
                    temp.s.push_back((coef*(radial1[j]-radial2[j])*coef2));          //scalar
                    temp.v.push_back((coef*(radial1[j]+radial2[j])*coef2));          //vector
                }
                append(temp,i);
                //cout<<temp.s[100]<<' '<<temp.v[100]<<endl;
            }
        }
       
        
    }
    
    
    void compute(int num){
        Solution x=solution[num];
        vector<int> kappas=x.kappas;
        for(int i=0;i<kappas.size();i++)
            for(int j=0;j<kappas.size();j++)
                compute_one_(i, j, num);
    }
    
    
//    void compute(int num){
//        Solution x=solution[num];    // The Solution class for 1 state
//        den temp;                    // a temp density for 1 channel
//        double m = x.m;              // m of this solution
//        double coef,coef2;
//        // coef is 3-j symbol, coef2 is just transformation between Legendre and Ylm
//        vector<int> kappas = x.kappas;        //all the kappa in the solution
//        int k1,k2;
//        vector<double> radial1,radial2;
//        vector<int> L_group;
//        for(int i=0;i<kappas.size();i++){
//            for(int j=0;j<kappas.size();j++){
//                //cout<<i<<' '<<j<<endl;
//                k1=x.wavefunctions[i].kappa;
//                k2=x.wavefunctions[j].kappa;
//                radial1=vector_multiple(x.wavefunctions[i].upper, x.wavefunctions[j].upper);
//                radial2=vector_multiple(x.wavefunctions[i].lower, x.wavefunctions[j].lower);
//                //cout<<x.wavefunctions[i].upper[100]<<' '<<x.wavefunctions[j].upper[100]<<radial1[100]<<endl;
//                for(int k=1;k<N;k++){
//                    radial1[k]/=pow(fx[k],2);
//                    radial2[k]/=pow(fx[k],2);
//                }
//                radial1[0]=radial1[1];
//                radial2[0]=radial2[1];
//                
//                L_group=get_possible_L(k1, k2);
//                for(int ii=0;ii<L_group.size();ii++){
//                    if (ii<max_L){
//                        coef=Angular_map.find(generate_key(m, k1, k2, ii))->second;
//                        //get the coef of angular term
//                        temp.s.clear();  //clear the vector
//                        temp.v.clear();
//                        for(int k=0;k<N;k++){
//                            coef2=sqrt((2.0*ii+1)/4/PI);
//                            temp.s.push_back((coef*radial1[k]-coef*radial2[k])/coef2);
//                            temp.v.push_back((coef*radial1[k]+coef*radial2[k])/coef2);
//                            //add the new value of this partial density to the density
//                            append(temp,ii);
//                        }
//                            
//                    }
//                }
//            }
//        }
//    }
    
    void Compute_all(){
        for(int i=0;i<solution.size();i++){
            compute(i);
        }
    }
};



//get density of everything for every iteration
void generate_density(vector<eig2> &occn, vector<eig2> &occp, vector<vector<double>> &dens,
                          vector<vector<double>> &denv,vector<vector<double>> &denp,vector<vector<double>> &den3){
    Density all_states_p(occp),all_states_n(occn);
    vector<den> density_p,density_n;
    vector<vector<double>> denvp,denvn,densn,densp;
    all_states_p.Compute_all();
    all_states_n.Compute_all();
    density_p=all_states_p.density;
    density_n=all_states_n.density;
    for(int i=0;i<max_L;i++){
        vector<double> temp1,temp2,temp3,temp4;
        for(int j=0;j<N;j++){
            temp1.push_back(density_p[i].s[j]);
            temp2.push_back(density_p[i].v[j]);
            temp3.push_back(density_n[i].s[j]);
            temp4.push_back(density_n[i].v[j]);
        }
        densp.push_back(temp1);
        denvp.push_back(temp2);
        densn.push_back(temp3);
        denvn.push_back(temp4);
    }
    
    for(int i=0;i<max_L;i++){
        for(int j=0;j<N;j++){
            dens[i][j]=densn[i][j]+densp[i][j];
            denv[i][j]=denvn[i][j]+denvp[i][j];
            den3[i][j]=denvp[i][j]-denvn[i][j];
            denp[i][j]=denvp[i][j];
        }
    }
}





#endif /* Density_h */
