#ifndef Bessel_h
#define Bessel_h

#include <iostream>
#include <math.h>
#include <vector>
#include "runge.h"
using namespace std;
//bessel functions with Imaginary
//x is a real number of r(mr< or mr>)
long double besselj(int l, double x){
    vector<long double> besseljs;
    if(l == 0){
        return sinh(x)/x;
    }
    if(l == 1){
        return cosh(x)/x - sinh(x)/pow(x,2);
    }
    besseljs.push_back(sinh(x)/x);                        //j0
    besseljs.push_back(cosh(x)/x - sinh(x)/pow(x,2));     //IM(j1)
    for(int i = 2; i <= l; i ++){
        //clate the ith bessel fucntions
        if (i% 2 == 0){
            long double a = (2*(i-1)+1.)*besseljs[i-1]/x - besseljs[i-2];
            besseljs.push_back(a);
        }
        else{
            // cout << i << endl;
            // cout << -(2*(i-1)+1.)*besseljs[i-1]/x << endl;
            // cout << -besseljs[i-2] << endl;
            long double a = -(2*(i-1)+1.)*besseljs[i-1]/x - besseljs[i-2];
            besseljs.push_back(a);
        }
    }
    return besseljs[l];
}

//hankel functions with imaginary
double hankel(int l, double x){
    vector<double> hankels;
    if(l==0){
        return -exp(-x)/x;
    }
    if(l==1){
        return exp(-x)/pow(x,2) + exp(-x)/x;
    }
    hankels.push_back(-exp(-x)/x);
    hankels.push_back(exp(-x)/pow(x,2) + exp(-x)/x);
    for(int i = 2; i <= l; i ++){
        if (i % 2 == 0){
            double a = (2*(i-1)+1.)*hankels[i-1]/x - hankels[i-2];
            hankels.push_back(a);
        }
        else{
            double a = -(2*(i-1)+1.)*hankels[i-1]/x - hankels[i-2];
            hankels.push_back(a);
        }
    }
    return hankels[l];
}

double htimesj(int l, double x){
    double a = besselj(l,x) * hankel(l,x);
    if( l%2 == 0){
        return a;
    }
    else{
        return -a;
    }
}

/*------------------------------------------------*/
//new implementations
double h_plus(int l, double r){
    int sign;
    if((l%2) == 0)
        sign = pow(-1,l+1+l/2);
    else
        sign = pow(-1,l+1+(l-1)/2);
    double a = exp(-r)/r;
    double b = 0;
    for(int i = 0; i <= l; i++){
        b += exp(lgamma(l+i+1) - lgamma(i+1) - lgamma(l-i+1) - i*log(2*r));
    }
    return sign*a*b;
}

double h_minus(int l, double r){
    double sign;
    if ((l%2) == 0)
        sign = pow(-1,l/2);
    else
        sign = pow(-1,(l-1)/2);
    double a = exp(r)/r;
    double b = 0;
    for(int i = 0; i <= l;i++){
        b += pow(-1,i)*exp(lgamma(l+i+1) - lgamma(i+1) - lgamma(l-i+1) - i*log(2*r));
    }
    return sign*a*b;
}

double besseljI(int l, double r){
    return (h_plus(l,r)+h_minus(l,r))/2.;
}


//returns j(l,Ix) * Ix
//for l == even, pure imaginary
//for l == odd,  pure real
double riccatijI(int l, double r){
    if((l%2)==0)
        return r*besseljI(l,r);
    else
        return -r*besseljI(l,r);
}


//return h(l,Ix) * Ix
//for l == even, pure imaginary
//for l == odd, pure real
double riccatihI(int l, double r){
    if((l%2)==0)
        return r*h_plus(l,r);
    else
        return -r*h_plus(l,r);
}

double friccatijI(double x, double y, double z, int l){
    return (1.0 + l*(l+1)/pow(x,2))*y;
}

double factorial2(int x){
    if(x==0){
        return 1;
    }
    if(x==1){
        return 1;
    }
    return (x)*factorial2(x-2);
}

vector<double> initbesselI(double x, int l){
    vector<double> ret(3, 0.0);
    if((l%2)==0){
        ret[0] = x;
        ret[1] = pow(-1,l/2)*pow(x,l+1)/factorial2(2*l+1);
        ret[2] = pow(-1,l/2)*(l+1)*pow(x,l)/factorial2(2*l+1);
    }
    else{
        ret[0] = x;
        ret[1] = pow(-1,(l+1)/2)*pow(x,l+1)/factorial2(2*l+1);
        ret[2] = pow(-1,(l+1)/2)*(l+1)*pow(x,l)/factorial2(2*l+1);
    }
    return ret;
}

vector<vector<double>> NormalizedRiccatijI(double h, int Numbers, int l){
    vector<double> initvector = initbesselI(h,l);
    vector<vector<double>> ret = runge_forward(initvector[0],initvector[1],initvector[2],h,Numbers,l,friccatijI);
    double norm = riccatijI(l, ret[0][Numbers])/ret[1][Numbers];
    for(int i = 0; i<= Numbers;i++){
        ret[1][i] = ret[1][i] * norm;
    }
    return ret;
}





#endif