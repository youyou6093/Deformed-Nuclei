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
#include <algorithm>
#include <unordered_map>
#include <string>
#include <iostream>
using namespace std;

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

int magic(int n){
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
    a.push_back(6);
    a.push_back(-5);
    a.push_back(4);
    a.push_back(-3);
    a.push_back(2);
    a.push_back(-1);
    a.push_back(-8);
    a.push_back(7);
    a.push_back(-6);
    a.push_back(5);
    a.push_back(-4);
    a.push_back(3);
    a.push_back(-2);
    a.push_back(1);
    int k = -1;
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
    return 2 * abs(min) - 1;
}

vector<string> generate_energy_level_map(){
    return {
        "0.-1", "0.-2", "0.1", "0.-3", "0.2", "1.-1", "0.-4", 
        "0.3", "1.-2", "1.1", "0.-5", "0.4", "1.-3", "1.2", 
        "2.-1", "0.-6"
    };
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


//map state 
unordered_map<string, vector<eig2>> generate_state_map(vector<eig2> &occ){
    unordered_map<string, vector<eig2>> ret;
    auto occ_object = get_solutions_object(occ);
    for (int i = 0; i < occ_object.size(); i++) {
        occ_object[i].get_primary_state();
        string key = to_string(occ_object[i].primary_state.n) + 
            "." + to_string(occ_object[i].primary_state.k); 
        if (!ret.count(key)) {
            vector<eig2> temp;
            temp.push_back(occ[i]);
            ret[key] = temp;
        } else {
            ret[key].push_back(occ[i]);
        }
     }
    // test
    // for (auto& item : ret) {
    //     cout << item.first << endl;
    //     auto temp = get_solutions_object(item.second);
    //     for (auto& state : temp) {
    //         state.get_primary_state();
    //         cout << state.primary_state << ";";
    //     }
    //     cout << endl;
    // }
    return ret;
} 


/* after we get all the raw solutions, change the solution to MEV and only keep the energies between -939 and 5,
 make sure all the energy is negative in the end*/
void get_solution(vector<eig2> &occp_raw,vector<eig2> &occn_raw,vector<eig2> &occp,vector<eig2> &occn){
    for(int i=0;i<occp_raw.size();i++){
        /*Transfer the energy units to Mev and also substract the mass*/
        occp_raw[i].solution.eigen_values = occp_raw[i].solution.eigen_values * hbarc - 939;
        if ((occp_raw[i].solution.eigen_values < 0) && (occp_raw[i].solution.eigen_values > -100)){
            occp.push_back(occp_raw[i]);
        }
    }
    for(int i=0;i<occn_raw.size();i++){
        occn_raw[i].solution.eigen_values=occn_raw[i].solution.eigen_values*hbarc-939;
        if ((occn_raw[i].solution.eigen_values < 0) && (occn_raw[i].solution.eigen_values>-100)){
            occn.push_back(occn_raw[i]);
        }
    }
    
    /*sort the solutions*/
    sort(occn.begin(),occn.end(),compare_eig2);
    sort(occp.begin(),occp.end(),compare_eig2);
    /*delete some solutions if I already got enough occ states*/
    if (neutron_number == 78) {
        vector<eig2> occn_temp;
        auto map = generate_state_map(occn);
        auto energy_level = generate_energy_level_map();
        for (int i = 0; i < energy_level.size(); i ++) {
            // cout << energy_level[i] << ' ' << map[energy_level[i]].size() << endl;
            for (auto state : map[energy_level[i]]) {
                occn_temp.push_back(state);
            }
        }
        occn = occn_temp;
    }
    if(occn.size()>neutron_number) occn.erase(occn.begin()+neutron_number,occn.end());
    if(occp.size()>proton_number) occp.erase(occp.begin()+proton_number,occp.end());
    
}



/* for specific m, return the solution eig2,contains m */
vector<eig2> get_temp_solution(vector<eig> &results,int m){
    eig2 temp;
    vector<eig2> temp_solution;
    for(int i=0;i<results.size();i++){
        temp.m=m;
        temp.solution=results[i];
        temp_solution.push_back(temp);
    }
    return temp_solution;
}



#endif /* utility_h */
