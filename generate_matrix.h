//
//  generate_matrix.h
//  deform_c++
//
//  Created by Junjie Yang on 4/18/17.
//  Copyright © 2017 Junjie Yang. All rights reserved.

#ifndef generate_matrix_h
#define generate_matrix_h
#include "simps.h"
#include "integrator.h"

extern int N;
extern int max_L;
extern vector<double> fx;
extern double hbarc,b,my_tolerance;
extern unordered_map<string, vector<vector<double>> > States_map;
extern unordered_map<string, double> Energy_map;
extern unordered_map<string, double > Angular_map;
extern const double PI;


/*Try to compute the angular dependece myself*/
extern"C" {
    double tj_(double* a1, double* a2, double* a3,double* b1, double* b2, double* b3);
    double sj_(double* j1, double* j2, double* j12,double* j3,double* j,double* j23);
    double coef9_(double* a1,double* a2,double* a3,
                  double* a4,double* a5,double* a6,
                  double* a7,double* a8,double* a9);
}

/*Compute the  cg coefficients based on the three-symbol*/
double cg_(double j1,double j2,double j3, double m1, double m2, double m3){
    return pow(-1, j1 - j2 - m3) * sqrt(2 * j3 + 1) * tj_(&j1, &j2, &j3, &m1, &m2, &m3);
}

/*Check the selection rule in alpha*/
int check(int l1, int l2, int L){
    int condition1,condition2,condition3;
    condition1 = (L >= abs(l1 - l2));
    condition2 = (L <= l1 + l2);
    condition3 = ((l1 + l2 + L) % 2 == 0);
    return condition1 && condition2 && condition3;
}

/* Alpha in jorge's notes
   These are 2j and 2m*/
double Angular_depedence(int j1, int j2, int l1, int l2, int m, int L){
    int sign = pow(-1, (m + 1) / 2);
    double a = sqrt(j1 + 1.0) * sqrt(j2 + 1.0);
    double b = cg_(j1/2.0, j2/2.0, L, m/2.0, -m/2.0, 0);
    double c = cg_(j1/2.0, j2/2.0, L, -0.5, 0.5, 0);
    return sign * a * b * c / 4.0 /PI;
}


/*Another Alpha but it uses kappa as input*/
double Angular_depedencek(int k1, int k2, int m, int L){
    int j1 = 2 * abs(k1) - 1;
    int j2 = 2 * abs(k2) - 1;
    int l1, l2;
    if (k1 > 0) l1 = k1;
    else l1 = -k1 - 1;
    if (k2 > 0) l2 = k2;
    else l2 = -k2 - 1;
    
    if (check(l1, l2, L) == 1)
        return Angular_depedence(j1, j2, l1, l2, m, L);
    else
        return 0;
}


/*Generat key for find the angular dependence
  right Now I didn't use it at all*/
string generate_key(double m,int k1,int k2,int L){
    if (m > 0)
        return to_string(int(m))+".5"+"."+to_string(k1)+"."+to_string(k2)+"."+to_string(L);
    else
        return "-" + to_string(int(abs(m)))+".5"+"."+to_string(k1)+"."+to_string(k2)+"."+to_string(L);
}

/* A struct that stores the params*/
struct params{
    int n1,n2,k1,k2,s1,s2;
    vector<double> g1,g2,f1,f2;
//    double *g1;
    int m;
    double diag_energy;
    vector<vector<double>> scalar_p,vector_p;
};




/*Calculate the matrix element for a specific L channel*/
double calculate_matrix_element(struct params & my_params,int L){
    double value1,value2,angular_term,value3;
    double E1,E2,E3;
//    vector<double> y1,y2,y3,y4,y5;
    double *y4 = new double[N];
//    vector<double> y4(N,0.0);
    /*The spherical symmetric potential part*/
    if (L==0){
        /*Diagonal part*/
        if((my_params.n1==my_params.n2)&&(my_params.s1==my_params.s2)&&(my_params.k1==my_params.k2))
            E1=my_params.diag_energy;
        else
            E1=0.0;
        /*Off diagnal part 1, only number of nodes affect this part*/
        if(my_params.k1==my_params.k2){
            for(int i=0;i<N;i++){
                /*hbarc here is to convert from Mev to fm-1*/
                value1=(my_params.scalar_p[L][i]+my_params.vector_p[L][i])*my_params.g1[i]*my_params.g2[i]/hbarc;
                value2=(-my_params.scalar_p[L][i]+my_params.vector_p[L][i])*my_params.f1[i]*my_params.f2[i]/hbarc;
                value3=-(1/pow(b,2))*(my_params.g1[i]*my_params.f2[i]+my_params.g2[i]*my_params.f1[i])*fx[i];
//                y3.push_back(value1+value2+value3);
                y4[i] = value1+value2+value3;
            }
            E2=simps(y4,&(fx[0]),N);
            // E2 = my_spline(y3, fx, my_tolerance).integral();
//            y3.clear();
        }
        else
            E2=0.0;
//        delete [] y4;
//        return E1+E2;
        E3 = E1 + E2;
    }
    
    /*The deformed part where kappa is not a good quantum number*/
    /*Be very careful here, the radial part potential here is P10
      Every potential here will be in terms of P, not Y */
    else{
        /* Angular part of the integral*/
        angular_term = 4 * PI / (2 * L + 1.0) * Angular_depedencek(my_params.k1, my_params.k2, my_params.m, L);
        if(angular_term == 0) E3=0.0;
        else{
            for(int i=0;i<N;i++){
                /* Radial part of the integral*/
                value1=angular_term * (my_params.scalar_p[L][i] + my_params.vector_p[L][i]) * my_params.g1[i] * my_params.g2[i] / hbarc;
                value2=angular_term * (-my_params.scalar_p[L][i] + my_params.vector_p[L][i]) * my_params.f1[i] * my_params.f2[i] / hbarc;
//                y4.push_back(value1+value2);
                y4[i] = value1+value2;
            }
            E3 = simps(y4,&(fx[0]),N);
            // E3 = my_spline(y4, fx, my_tolerance).integral();
//            y4.clear();
        }
    }
    delete [] y4;
    return E3;
}



/*this function,get M[i,j] where M is the final matrix*/
void generate_matrix(vector<vector<double>> &M,vector<vector<double>> &scalar_potential,vector<vector<double>> & vector_potential,vector<State> &states,int m,int i,int j){
    unordered_map<string, vector<vector<double>>>::iterator States_ptr;
    unordered_map<string, double>::iterator Energy_ptr;
    
    double matrix_element=0.0;
    struct params my_params;
    my_params.m=m;
    my_params.n1=states[i].n;
    my_params.n2=states[j].n;
    my_params.k1=states[i].k;
    my_params.k2=states[j].k;
    my_params.s1=states[i].sign;
    my_params.s2=states[j].sign;
    my_params.diag_energy=Energy_map.find(states[i].key)->second; //the energy for states[i]
    
    
    States_ptr=States_map.find(states[i].key);
//    my_params.g1=&(States_ptr->second[0][0]);
    my_params.g1=States_ptr->second[0];
    my_params.f1=States_ptr->second[1];
    States_ptr=States_map.find(states[j].key);
    my_params.g2=States_ptr->second[0];
    my_params.f2=States_ptr->second[1];
    my_params.scalar_p=scalar_potential;
    my_params.vector_p=vector_potential;
    

   
    for(int L=0; L <max_L; L++){
        matrix_element+=calculate_matrix_element(my_params, L);
    }
    
    M[i][j]=matrix_element;
}
    

/*generate the full matrix for a specific m and a specific particle type
 (this is determined by the scalar or vector potential....*/
vector<vector<double>> generate_full_matrix(vector<vector<double>> &scalar_potential,vector<vector<double>> & vector_potential,vector<State> &states,int m){
    int matrix_size=int(states.size());
    vector<vector<double>> M(matrix_size,vector<double>(matrix_size,0.0));
#pragma omp parallel for
    for(int i=0;i<states.size();i++){
        for(int j=i;j<states.size();j++){   //this is half of the matrix
            generate_matrix(M, scalar_potential, vector_potential, states, m, i, j);
        }
    }
    for(int j=0;j<states.size();j++)
        for(int i=0;i<j;i++){
            M[j][i]=M[i][j];               //the other half has already been computed
        }
    return M;
}

#endif /* generate_matrix_h */
