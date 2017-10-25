//
//  effective_density.h
//  deform_c++
//
//  Created by Junjie Yang on 7/17/17.
//  Copyright © 2017 Junjie Yang. All rights reserved.
//


#ifndef effective_density_h
#define effective_density_h
#include "Green-method.h"


extern"C" {
    double tj_(double* a1, double* a2, double* a3,double* b1, double* b2, double* b3);
    double sj_(double* j1, double* j2, double* j12,double* j3,double* j,double* j23);
    double coef9_(double* a1,double* a2,double* a3,
                  double* a4,double* a5,double* a6,
                  double* a7,double* a8,double* a9);
}         //declare that there is something in the .o library




extern vector<double> empty;
extern double lambdas,lambdav,lambda,ka,ks,ms,mv,mp,mg,gs,gv,gp,gg;



double pcg_(double j1,double j2,double j3){            //compute the parity cg coefficients based on the three-symbol
    double m1=0;
    double m2=0;
    double m3=0;
    return pow(-1,j1-j2)*sqrt(2*j3+1)*tj_(&j1, &j2, &j3, &m1,&m2,&m3);
}
//need to verify


vector<int> possible_L(int l1,int l2){
    vector<int> Ls;
    for(int i=abs(l1-l2);i<l1+l2+1;i++){
        if ((i+l1+l2) % 2 == 0){
            Ls.push_back(i);
        }
    }
    return Ls;
}

vector<vector<double>> multiplication(vector<vector<double>> &vec1,vector<vector<double>> &vec2){
    vector<vector<double>> new_vec;
    vector<double> temp;
    vector<int> all_L;
    int Big_L;
    for(int i=0;i<max_L;i++){
        new_vec.push_back(empty);
    }
    for(int i=0;i<max_L;i++){
        for(int j=0;j<max_L;j++){
            all_L=possible_L(i, j);             //get all possible  L with i and j
            for(int k=0;k<all_L.size();k++){    //loop through all possible L
                Big_L=all_L[k];
                if (Big_L<max_L){               //if valid
                    temp.clear();
                    for(int ii=0;ii<N;ii++){    //get temp
                        temp.push_back(vec1[i][ii]*vec2[j][ii]*pow(pcg_(i, j, Big_L),2));
                    }
                    for(int ii=0;ii<N;ii++){
                        new_vec[Big_L][ii]+=temp[ii];  //update the new vector
                    }
                }
            }
        }
    }
    
    return new_vec;
}


//get the effective density
void get_effective_density(vector<vector<double>> &EFF_Phi,vector<vector<double>>  &EFF_B,vector<vector<double>> &EFF_A,vector<vector<double>> &EFF_W,
                           vector<vector<double>> &Phi,vector<vector<double>> &W,vector<vector<double>> &B,vector<vector<double>> &A,
                           vector<vector<double>> &dens,vector<vector<double>> &denv,vector<vector<double>> & den3,vector<vector<double>> denp){
    
    vector<vector<double>> W2,W3,Phi2,Phi3,B2,po1,po2,po3,po4;
    W2=multiplication(W, W);
    W3=multiplication(W2, W);
    Phi2=multiplication(Phi, Phi);
    Phi3=multiplication(Phi2, Phi);
    B2=multiplication(B, B);
    po1=multiplication(B2, Phi);
    po2=multiplication(B2, W);
    po3=multiplication(Phi2, B);
    po4=multiplication(W2, B);
    
    for(int i=0;i<max_L;i++){
        for(int j=0;j<N;j++){
            EFF_Phi[i][j] = dens[i][j] + 2*lambdas*po1[i][j]/pow(hbarc,3)-(lambda/6.0)*Phi3[i][j]/pow(hbarc,3)-(ka/2.0)*Phi2[i][j]/pow(hbarc,2);
            EFF_W[i][j] = denv[i][j] - 2*lambdav*po2[i][j]/pow(hbarc,3) - (ks/6.0)*W3[i][j]/pow(hbarc,3);
            EFF_B[i][j] = 0.5*(den3[i][j] - 4*lambdas*po3[i][j]/pow(hbarc,3) - 4*lambdav*po4[i][j]/pow(hbarc,3));
        }
    }
    
}


void get_potential(vector<vector<double>> &EFF_Phi,vector<vector<double>>  &EFF_B,vector<vector<double>> &EFF_A,vector<vector<double>> &EFF_W,
                   vector<vector<double>> &Phi,vector<vector<double>> &W,vector<vector<double>> &B,vector<vector<double>> &A,
                   vector<vector<double>> &dens,vector<vector<double>> &denv,vector<vector<double>> & den3,vector<vector<double>> denp,
                   int L){
    for(int i=1;i<N;i++){
        Phi[L][i] = 2 * Phi[L][i]/3.0 + hbarc * gs * klein(ms,  EFF_Phi[L] , i)/3.0; //not sure whether I can just put gs outside
        W[L][i] = 2 * W[L][i]/3.0 + hbarc * gv * klein(mv, EFF_W[L], i)/3.0;
        B[L][i] = 2 * B[L][i]/3.0 + hbarc * gp * klein(mp , EFF_B[L], i)/3.0;
        A[L][i] = hbarc * gg * klein(mg, denp[L], i);
    }
    Phi[L][0] = Phi[L][1];
    W[L][0] = W[L][1];
    B[L][0] = B[L][1];
    A[L][0] = A[L][1];
    
    
    
}

void update_potential(vector<vector<double>> &EFF_Phi,vector<vector<double>>  &EFF_B,vector<vector<double>> &EFF_A,vector<vector<double>> &EFF_W,vector<vector<double>> &Phi,vector<vector<double>> &W,vector<vector<double>> &B,vector<vector<double>> &A,vector<vector<double>> &dens,vector<vector<double>> &denv,vector<vector<double>> & den3,vector<vector<double>> denp){
    for(int i=0;i<5;i++){        //how many iteration I want to solve klein-gordon equation
        get_effective_density(EFF_Phi, EFF_B, EFF_A, EFF_W, Phi, W, B, A, dens, denv, den3, denp);
        for(int j=0;j<max_L;j++){
            get_potential(EFF_Phi, EFF_B, EFF_A, EFF_W, Phi, W, B, A, dens, denv, den3, denp, j);
        }
    }
}


#endif /* effective_density_h */
