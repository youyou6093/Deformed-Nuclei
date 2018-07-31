#ifndef Vector_compute_h
#define Vector_compute_h
//tested

#include<vector>
using namespace std;

vector<double> vector_add(vector<double> & a,vector<double>& b){
    vector<double> c;
    for(int i=0;i<a.size();i++){
        c.push_back(a[i]+b[i]);
    }
    return c;
}

vector<double> vector_multiple(vector<double> &a,vector<double> &b){
    vector<double> c(a.size(), 0.0);
    for(int i=0;i<a.size();i++) c[i] = a[i] * b[i];
    return c;
}

vector<double> vector_time(vector<double> &a, double & b){
    vector<double> c;
    for(int i=0;i<a.size();i++){
        c.push_back(a[i]*b);
    }
    return c;
}




#endif
