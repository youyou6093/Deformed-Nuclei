#include <iostream>
#include <math>
//#include "generate_matrix.h"
using namespace std;


extern"C" {
    double tj_(double* a1, double* a2, double* a3,double* b1, double* b2, double* b3);
    double sj_(double* j1, double* j2, double* j12,double* j3,double* j,double* j23);
    double coef9_(double* a1,double* a2,double* a3,
                  double* a4,double* a5,double* a6,
                  double* a7,double* a8,double* a9);
}         //declare that there is something in the .o library

//compute the  cg coefficients based on the three-symbol
double cg_(double j1,double j2,double j3, double m1, double m2, double m3){
    return pow(-1, j1 - j2 - m3) * sqrt(2 * j3 + 1) * tj_(&j1, &j2, &j3, &m1, &m2, &m3);
}


int check(int l1, int l2){
    int condition1,condition2,condition3;
    condition1 = (1 >= abs(l1 - l2));
    condition2 = (1 <= l1 + l2);
    condition3 = ((l1 + l2 + 1) % 2 == 0);
    return condition1 && condition2 && condition3;
}

/* alpha in jorge's notes*/
/* these are 2j and 2m*/
double Angular_depedence(int j1, int j2, int l1, int l2, int m, int L){
    int sign = pow(-1, (m + 1) / 2);
    double a = sqrt(j1 + 1.0) * sqrt(j2 + 1.0);
    double b = cg_(j1/2.0, j2/2.0, L, m/2.0, -m/2.0, 0);
    double c = cg_(j1/2.0, j2/2.0, L, -0.5, 0.5, 0);
    return sign * a * b * c / 4.0 /PI;
}

double Angular_depedencek(int k1, int k2, int m, int L){
    int j1 = 2 * abs(k1) - 1;
    int j2 = 2 * abs(k2) - 1;
    int l1, l2;
    if (k1 > 0) l1 = k1;
    else l1 = -k1 - 1;
    if (k2 > 0) l2 = k2;
    else l2 = -k2 - 1;
    return Angular_depedence(j1, j2, l1, l2, m, L);
}








int main(){
    cout << Angular_depedencek(-14,-14,-27,0) << endl;
}
