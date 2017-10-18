//  diagnolize a normal matrix
//  symm.h
//  gsl_matrix
//
//  Created by Junjie Yang on 5/9/17.
//  Copyright © 2017 Junjie Yang. All rights reserved.
// tested

#ifndef symm_h
#define symm_h

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include<vector>
#include<iostream>
using namespace std;
//be careful, this program will change the content of the array passed in if you use the array constructor

struct eig{              //a structure that stores the eigenvectors
    double eigen_values;
    vector<double> eigen_vectors;
    
};

bool compare_eig(eig a,eig b){
    return (a.eigen_values<b.eigen_values);
}


class matrix_diag{
private:
    double *m_array;             //1d array
    //--------------------------------------
    gsl_matrix_view m;             //gsl matrix
    gsl_vector *eval = NULL;              //eigenvalues later
    gsl_matrix *evec = NULL;              //eigenvectors later
    gsl_eigen_symmv_workspace * w;
    //-----------------------------------------
    
    void diagonalize(){
        m=gsl_matrix_view_array(m_array, size, size);
        eval=gsl_vector_alloc(size);
        evec=gsl_matrix_alloc(size, size);
        w=gsl_eigen_symmv_alloc(size);
        gsl_eigen_symmv (&m.matrix, eval, evec, w);
        gsl_eigen_symmv_free(w);
        gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);
    }

    
public:
    
    vector<eig> results;         //eigenvalues eigenvectors.
    double size;                 //size of the matrix(not total size, just the dimension of the matrix)
    
    ~matrix_diag(){              //destructor
        gsl_vector_free(eval);
        gsl_matrix_free(evec);
    }
    
    void get_results(){             //get eigenvalues-eigenvectors pair(real part)
        diagonalize();              //diagnolize the matrix
        eig temp;                   //stores a pair of eigenvalues and eigenvectors
        for(int i=0;i<size;i++){
            double eval_i=gsl_vector_get(eval,i);                   //get ith eigenvalues
            gsl_vector_view evec_i=gsl_matrix_column(evec,i);       //get ith eigenvectors
            temp.eigen_values=eval_i;                               //store eigenvalues
            for(int j=0;j<size;j++){
                double x = gsl_vector_get(&evec_i.vector,j);
                temp.eigen_vectors.push_back(x);                   //store ith eigenvectors,
                                                                   //jth values
            }
            results.push_back(temp);
            temp.eigen_vectors.clear();                            //clear the temp;
        }
    }
    
    
    matrix_diag(){
        m_array=0;
        size=0;
    }
    
    matrix_diag(double *data,double n){    //initialize the matrix from array
        //pass by reference,will change data itself
        m_array=data;
        size=n;
    }
    
    
    matrix_diag(vector<double> &data,double n){           //it is very important to pass by reference!!!
        m_array=&(data[0]);
        size=n;
    }
    
    
    //this one is not working, I'll change that later
    matrix_diag(vector<vector<double>> &M,int size){   //initialize the matrix from 2d vector
        
        
        
        //won't change M itself,deep copy M to m_array
        //cout<<M.size()<<' '<<size<<endl;
        double a[size*size];
        int ptr=0;
        for(int i=0;i<size;i++){
            for(int j=0;j<size;j++){
                a[ptr]=M[i][j];
                ptr++;
            }
        }
        cout<<ptr<<' '<<size*size<<endl;
        this->size=size;
        m_array=a;
    }
    
    
    
};



#endif /* symm_h */