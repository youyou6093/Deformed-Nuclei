//  diagnolize a normal matrix
//  symm.h
//  Created by Junjie Yang on 5/9/17.
//  Copyright © 2017 Junjie Yang. All rights reserved.

#ifndef symm_h
#define symm_h

#include <stdio.h>
#include<vector>
#include<iostream>
#include "Eigen/Eigenvalues"
using namespace Eigen;
using namespace std;

struct eig{              //a structure that stores the eigenvectors
    double eigen_values;
    vector<double> eigen_vectors;
    
};

bool compare_eig(eig a,eig b){
    return (a.eigen_values<b.eigen_values);
}


class matrix_diag{
private:
    MatrixXd m;
//    VectorXd eval;
//    MatrixXd evec;
    SelfAdjointEigenSolver<MatrixXd> eigenSolver;
    void diagonalize(){
        eigenSolver.compute(m);
    }    
public:
    
    vector<eig> results;         //eigenvalues eigenvectors.
    int size;                    //size of the matrix(not total size, just the dimension of the matrix)
    
    ~matrix_diag(){              //destructor
    }
    
    void get_results(){             //get eigenvalues-eigenvectors pair(real part)
        results.clear();
        diagonalize();              //diagnolize the matrix
        eig temp;                   //stores a pair of eigenvalues and eigenvectors
        for(int i = 0; i < size; i++){
            temp.eigen_values = eigenSolver.eigenvalues()[i];   //ith eigenvalue
            for(int j =0;j<size;j++){
                temp.eigen_vectors.push_back(eigenSolver.eigenvectors().col(i)[j]);
            }   //the jth element of ith eigenvalues
            results.push_back(temp);
            temp.eigen_vectors.clear();
        }
    }
    
    
    matrix_diag(vector<vector<double>> &M){  
        size = M.size();
        m = MatrixXd::Zero(size,size);
        for (size_t i = 0; i < size; ++i)
            for (size_t j = 0; j < size; ++j)
                m(i,j) = M[i][j];
    }
};

#endif /* symm_h */
