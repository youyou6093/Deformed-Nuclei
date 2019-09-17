#include <iostream>
#include <vector>
#include <fstream>
using namespace std;
double hbarc = 197.326;
double ms=491.500/hbarc;
double mv=782.5/hbarc;
double mp=763.0/hbarc;
double b=2.4,rmin=0.0,rmax=20.0;
int N=801;
double h=(rmax-rmin)/(N-1);       //grid size
vector<double> fx;                //x axis

#include "Bessel.h"
#include "runge.h"
#include "integrator.h"

void Get_bessels(int type, my_spline & riccatijIs, vector<vector<double>> &JIs, vector<vector<double>> &HIs, int L){
    double m;
    if(type == 0){
        m = ms;
    }
    else if(type == 1){
        m = mv;
    }
    else if(type == 2){
        m = mp;
    }
    for(int i = 1; i < N; i++){
        if ((fx[i] * m) < 1){
            JIs[type][i] = riccatijIs.eval(fx[i] * m);
        }
        else{
            JIs[type][i] = riccatijI(L, fx[i] * m);
        }
        HIs[type][i] = riccatihI(L, fx[i] * m);
    }
    JIs[type][0] = 0.0;
    HIs[type][0] = 0.0;
}


int main(){
    int L;
    cin >> L;
    for(int i = 0; i < N; i++) fx.push_back(h*i + rmin);
    vector<vector<double>> JIs(3, vector<double>(N, 0.0)), HIs(3, vector<double>(N, 0.0));
    double h_bessel = 0.001;
    int Numbers = 1000;
    vector<vector<double>> ret = NormalizedRiccatijI(h_bessel, Numbers, L);
    vector<double> ydata(Numbers + 1, 0.0);
    vector<double> xdata(Numbers + 1, 0.0);
    for(int i = 0; i <= Numbers; i++){
        ydata[i] = ret[1][i];
        xdata[i] = ret[0][i];
    }
    my_spline riccatijIs = my_spline(ydata, xdata, 0.001);
    Get_bessels(0, riccatijIs, JIs, HIs, L);
    Get_bessels(1, riccatijIs, JIs, HIs, L);
    Get_bessels(2, riccatijIs, JIs, HIs, L);
    ofstream os1,os2;
    os1.open("bessel.dat");
    os2.open("hankel.dat");
    for(int i = 0; i < N; i++){
        os1 << fx[i] << ' ' << JIs[0][i] << ' ' << JIs[1][i] << ' ' << JIs[2][i] << endl;
        os2 << fx[i] << ' ' << HIs[0][i] << ' ' << HIs[1][i] << ' ' << HIs[2][i] << endl;
    }
    os1.close();
    os2.close();
    
    
    return 0;
}
