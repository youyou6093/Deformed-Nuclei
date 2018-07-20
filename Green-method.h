//
//  Green-method.h
//  deform_c++
//
//  Created by Junjie Yang on 7/17/17.
//  Copyright © 2017 Junjie Yang. All rights reserved.


#ifndef Green_method_h
#define Green_method_h
#include "integrator.h"
#include "Bessel.h"
#include <thread>
#include <chrono>
using namespace std;

double klein2(double mass, vector<double> density, int index, int L, my_spline & riccatijIs, vector<vector<double>> JIs, vector<vector<double>> HIs, int type){
    double r;
    r = fx[index]; //the point that I want to know the potential
    double coef = 1/(mass*r);
    double coef1,coef2;
    vector<double> inte1, inte2, fx1, fx2;
    coef1 = riccatihI(L, r * mass);
    if(r * mass < 1){
        coef2 = riccatijIs.eval(r*mass);
    }
    else{
        coef2 = riccatijI(L, r*mass);
    }

    for(int i = 0; i < index + 1;i ++){  //left half
//        double a;
//        if (i == 0){
//            a = 0;
//            // b = 0;
//        }
//        else if((fx[i]*mass)< 1){
//            a = riccatijIs.eval(fx[i]*mass);       //evaluate j^(imx')
//            // b = riccatihI(L, r*mass);               //evaluate h^(imx)
//        }
//        else{
//            a = riccatijI(L, fx[i]*mass);
//            // b = riccatihI(L, r*mass);
//        }
        
        
        if((L % 2) == 0)   //need an extra minus because the even channels are imaginary
            inte1.push_back(-JIs[type][i]*coef1*fx[i]*density[i]); //the left integral
        else
            inte1.push_back(JIs[type][i]*coef1*fx[i]*density[i]);
        fx1.push_back(fx[i]);
    }
    for(int i = index; i < N; i++){
        double a;
//        if(r*mass < 1){
//            // a = riccatijIs.eval(r*mass);         //evaluate j^(imx)
//            a = riccatihI(L, fx[i]*mass);         //evaluate h^(imx')
//        }
//        else{
//            // a = riccatijI(L, r*mass);
//            a = riccatihI(L, fx[i]*mass);
//        }

        if((L % 2) == 0)
            inte2.push_back(-coef2*HIs[type][i]*fx[i]*density[i]);
        else
            inte2.push_back(coef2*HIs[type][i]*fx[i]*density[i]);
        fx2.push_back(fx[i]);
    }


    //manual spline for extreme case
    if(fx1.size()==2){
        double manual_spline_y = 0.5*(inte1[0] + inte1[1]);
        double manual_spline_x = 0.5*(fx1[0] + fx1[1]);
        fx1.insert(fx1.begin()+1, manual_spline_x);
        inte1.insert(inte1.begin()+1, manual_spline_y);
    }
    if(fx2.size()==2){
        double manual_spline_y = 0.5*(inte2[0]+inte2[1]);
        double manual_spline_x = 0.5*(fx2[0]+fx2[1]);
        fx2.insert(fx2.begin()+1,manual_spline_x);
        inte2.insert(inte2.begin()+1,manual_spline_y);
    }

    // cout << mass << ' '  << index << ' ' << fx2[0] << ' ' << inte2[0]*coef << ' ' << density[index] << endl;

    // for(int i = 0; i < fx[1].size; i++){
    // 	cout << 
    // }


    // std::this_thread::sleep_for(std::chrono::milliseconds(2000));

    if(fx2.size() == 1){
        return coef*my_spline(inte1,fx1,my_tolerance).integral();
    }
    else{
        return coef*(my_spline(inte1,fx1,my_tolerance).integral() + my_spline(inte2,fx2,my_tolerance).integral());
    }



}

/*need more checking*/
double poisson(vector<double> density, int index, int L){
    double r = fx[index];
    double coef = 1/(2*L+1.);
    vector<double> inte1,inte2,fx1,fx2;
    for(int i = 0; i < index + 1; i++){
        double a = pow(fx[i],2) * density[i];
        inte1.push_back(a*pow(fx[i],L)/pow(r,L+1));
        fx1.push_back(fx[i]);
    }
    for(int i = index; i < N;i++){
        double a = pow(fx[i],2)*density[i];
        inte2.push_back(a*pow(r,L)/pow(fx[i],L+1));
        fx2.push_back(fx[i]);
    }
    
    if(fx1.size()==2){           //in this situation, we will be not able to perform a  cubic spline
        double manual_spline_y = 0.5 * (inte1[0]+inte1[1]);
        double manual_spline_x = 0.5 * (fx1[0] + fx1[1]);
        //cout<<manual_spline_x<<endl;
        fx1.insert(fx1.begin()+1,manual_spline_x);
        inte1.insert(inte1.begin()+1, manual_spline_y);
    }
    if(fx2.size()==2){          //in this situation, we will be not able to perform a  cubic spline
        double manual_spline_y = 0.5*(inte2[0]+inte2[1]);
        double manual_spline_x = 0.5*(fx2[0]+fx2[1]);
        fx2.insert(fx2.begin()+1, manual_spline_x);
        inte2.insert(inte2.begin()+1, manual_spline_y);
    }
    
    if(fx2.size()==1)   //the second part of the integral is essentially zero
        return coef*my_spline(inte1,fx1,my_tolerance).integral();
    else
        return coef*(my_spline(inte1,fx1,my_tolerance).integral()+my_spline(inte2,fx2,my_tolerance).integral());
    
}

double klein(double mass,vector<double> density,int index){
    double r,coef1,coef2,coef3,temp1,temp2;
    vector<double> inte1,inte2,fx1,fx2;
    r=fx[index];
    coef1=1.0/(2*mass);
    coef2=exp(-mass*r)/r;
    coef3=(exp(mass*r)-exp(-mass*r))/r;
    for(int i=0;i<index+1;i++){
        temp1=fx[i]*(exp(mass*fx[i])-exp(-mass*fx[i]))*density[i];
        inte1.push_back(temp1);
        fx1.push_back(fx[i]);
    }
    for(int i=index;i<N;i++){
        temp2=fx[i]*exp(-mass*fx[i])*density[i];
        inte2.push_back(temp2);
        fx2.push_back(fx[i]);
    }
    if(fx1.size()==2){           //in this situation, we will be not able to perform a  cubic spline
        double manual_spline_y = 0.5 * (inte1[0]+inte1[1]);
        double manual_spline_x = 0.5 * (fx1[0] + fx1[1]);
        //cout<<manual_spline_x<<endl;
        fx1.insert(fx1.begin()+1,manual_spline_x);
        inte1.insert(inte1.begin()+1, manual_spline_y);
    }
    if(fx2.size()==2){          //in this situation, we will be not able to perform a  cubic spline
        double manual_spline_y = 0.5*(inte2[0]+inte2[1]);
        double manual_spline_x = 0.5*(fx2[0]+fx2[1]);
        fx2.insert(fx2.begin()+1, manual_spline_x);
        inte2.insert(inte2.begin()+1, manual_spline_y);
    }

    // cout << mass << ' '  << index << ' ' << fx2[0] << ' ' << coef1*coef3*inte2[0] << ' ' << density[index] << endl;
    // std::this_thread::sleep_for(std::chrono::milliseconds(2000));

    if(fx2.size()==1)   //the second part of the integral is essentially zero
        return coef1*(coef2*my_spline(inte1,fx1,my_tolerance).integral());
    else
        return coef1*(coef2*my_spline(inte1,fx1,my_tolerance).integral()+coef3*my_spline(inte2,fx2,my_tolerance).integral());
    
}



#endif /* Green_method_h */
