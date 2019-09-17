#ifndef Runge_h
#define Runge_h

#include<iostream>
#include<vector>
using namespace std;

vector<vector<double>> runge_forward(double x0,double y0,double z0,double h,int Numbers,int l,double (*f)(double, double, double, int)){
    vector<vector<double>> ret(3, vector<double>(Numbers + 1, 0.0));
    ret[0][0] = x0;
    ret[1][0] = y0;
    ret[2][0] = z0;
    double x = x0;
    double y = y0;
    double z = z0;
    double k11,k12,k21,k22,k31,k32,k41,k42;
    for(int i = 0; i < Numbers; i ++){
        k11 = z;
        k12 = (*f)(x,y,z,l);
        k21 = z+(h/2)*k12;
        k22 = (*f)(x+h/2,y+(h/2)*k11,z+(h/2)*k12,l);
        k31 = z+(h/2)*k22;
        k32 = (*f)(x+h/2,y+(h/2)*k21,z+(h/2)*k22,l);
        k41 = z+h*k32;
        k42 = (*f)(x+h,y+h*k31,z+h*k32,l);
        y = y+(h/6)*(k11+2*k21+2*k31+k41);
        z = z+(h/6)*(k12+2*k22+2*k32+k42);
        x = x+h;
        ret[0][i+1] = x;
        ret[1][i+1] = y;
        ret[2][i+1] = z;
    }
    return ret;
}

#endif
