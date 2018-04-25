//tested
//  simps.h
//  deform_c++
//
//  Created by Junjie Yang on 5/4/17.
//  Copyright © 2017 Junjie Yang. All rights reserved.
//

#ifndef simps_h
#define simps_h

#include<iostream>
#include<math.h>
#include<vector>
using namespace std;


double simps(double *y,double *x, int size){
//    int size=int(x.size());
    if (size%2==0){
        cout<<"size has to be odd"<<endl;
        return -1.0;
    }
    double a=x[0],b=x[size-1],h=(b-a)/(size-1);
    double sum=0.0;
    for(int i=1;i<size-1;i+=2) sum+=4.0*y[i];
    for(int i=2;i<size-2;i+=2) sum+=2.0*y[i];
    return h/3.0*(y[0]+y[size-1]+sum);
}




double simps(vector<double> &y,vector<double> &x){
    int size=int(x.size());
    if (size%2==0){
        cout<<"size has to be odd"<<endl;
        return -1.0;
    }
    double a=x[0],b=x[size-1],h=(b-a)/(size-1);
    double sum=0.0;
    for(int i=1;i<size-1;i+=2) sum+=4.0*y[i];
    for(int i=2;i<size-2;i+=2) sum+=2.0*y[i];
    return h/3.0*(y[0]+y[size-1]+sum);
}




#endif /* simps_h */
