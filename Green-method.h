//
//  Green-method.h
//  deform_c++
//
//  Created by Junjie Yang on 7/17/17.
//  Copyright © 2017 Junjie Yang. All rights reserved.


#ifndef Green_method_h
#define Green_method_h
#include "integrator.h"



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
    if(fx2.size()==1)   //the second part of the integral is essentially zero
        return coef1*(coef2*my_spline(inte1,fx1,my_tolerance).integral());
    else
        return coef1*(coef2*my_spline(inte1,fx1,my_tolerance).integral()+coef3*my_spline(inte2,fx2,my_tolerance).integral());
    
}



#endif /* Green_method_h */
