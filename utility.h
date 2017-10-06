//
//  utility.h
//  deform_c++
//
//  Created by Junjie Yang on 8/3/17.
//  Copyright © 2017 Junjie Yang. All rights reserved.
//

#ifndef utility_h
#define utility_h
#include "integrator.h"

extern int proton_number;
extern int neutron_number;
extern double hbarc,ms,mv,mp,mg,gs,gv,gp,gg,lambdas,lambdav, lambda, ks,ka;
double compute_energy(vector<eig2> &occp,vector<eig2> &occn,vector<vector<double>> &Phi,
                      vector<vector<double>> &W, vector<vector<double>> &B, vector<vector<double>> &A,
                      vector<vector<double>> &dens, vector<vector<double>> &denv, vector<vector<double>> &den3,
                      vector<vector<double>> &denp){
    double ecm=3.0*41.0*pow((proton_number+neutron_number),-4.0/3)/4.0;
    double total_energy_part1=0.0;
    for(int i=0;i<occn.size();i++)
        total_energy_part1+=occn[i].solution.eigen_values;
    for(int i=0;i<occp.size();i++)
        total_energy_part1+=occp[i].solution.eigen_values;
    //vector<double> inte1,inte2,inte3,inte4,part1,part2,part3,part4;
    vector<double> inte;
    for(int i=0;i<N;i++){
        double temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,Final;
        temp1 = - Phi[0][i] * dens[0][i];
        temp2 = W[0][i] * denv[0][i];
        temp3 = 0.5 * B[0][i] * den3[0][i];
        temp4 = A[0][i] * denp[0][i];
        temp5 = 2 * PI * pow(fx[i],2) *(temp1+temp2+temp3+temp4);
        temp6 = 1.0/3 * PI * pow(fx[i],2) *(ka*hbarc*pow(Phi[0][i],3)+0.5*lambda*pow(Phi[0][i],4))/pow(hbarc,3);
        temp7 = -4 * PI * pow(fx[i],2) * ks* pow(W[0][i],4)/24.0/pow(hbarc,3);
        temp8 = -4*PI*pow(fx[i],2)*(lambdas*pow(Phi[0][i],2)+lambdav*pow(W[0][i],2)*pow(B[0][i],2))/pow(hbarc,3);
        Final = temp5+temp6+temp7+temp8;
        inte.push_back(Final);
    }
    return (total_energy_part1 - my_spline(inte,fx,my_tolerance).integral())/(proton_number+neutron_number) - ecm;
    
    
    
}

double magic(int n){
    vector<int> a;
    a.push_back(-1);
    a.push_back(-2);
    a.push_back(1);
    a.push_back(-3);
    a.push_back(-1);
    a.push_back(2);
    a.push_back(-4);
    a.push_back(-2);
    a.push_back(3);
    a.push_back(1);
    a.push_back(-5);
    a.push_back(4);
    a.push_back(-3);
    a.push_back(2);
    a.push_back(-1);
    a.push_back(-6);
    a.push_back(5);
    a.push_back(-4);
    a.push_back(3);
    a.push_back(-2);
    a.push_back(+1);
    a.push_back(-7);
    int k;
    for(int i=0;i<30;i++){
        n-=2*abs(a[i]);
        if(n<=0){
            k=i;
            break;
        }
    }
    int min=-1;
    for(int i=0;i<k+1;i++){
        if (a[i]<min){
            min=a[i];
        }
    }
    return abs(min)-0.5;
}






#endif /* utility_h */
