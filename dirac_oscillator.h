//tested

#include "HOStuff.h"
#include<iostream>
using namespace std;
//becareful this file returns radial part

class dirac_oscillator{
private:
    int p_sign;
    int p_n;
    int p_kappa;
    double p_b,p_h;
    double p_rmin;
    double p_rmax;
    double p_N;
    double M;
    void transform(int n,int l,double nu,double norm,vector<double> & wave){
        HOLevelWFs ho(n,l,nu,p_rmin,p_rmax,p_N);
        wave=ho.R;
        for(int i=0;i<p_N;i++)
            wave[i]=norm*wave[i]*(p_rmin+p_h*i);
    }
    
public:
    double E;
    vector<double> upper;
    vector<double> lower;
    
    
    
    dirac_oscillator(int sign,int n,int kappa,double b,double rmin,double rmax,int N){
        //the parameters of the ho_basis
        M=939.0/197.326;
        p_sign=sign;
        p_n=n;
        p_kappa=kappa;
        p_b=b;
        p_rmin=rmin;
        p_rmax=rmax;
        p_N=N;
        // finish this part
        p_h=(rmax-rmin)/(N-1);
        int n2,l,l2;
        double alpha,norm;
        double nu=0.5/pow(b,2);
        //cout<<sign<<endl;
        
        if (sign>0){       //positive energy states
            //cout<<"big"<<endl;
            if (kappa<0){
                E=pow((4*n/pow(b,2)+pow(M,2)),0.5);
                n2=n-1;
                l=-kappa-1;
                l2=l+1;
                alpha=-2.0*pow(n,0.5)/(b*(E+M));
                norm=pow((1.0/(1+pow(alpha,2))),0.5);   //normalizations
                
                transform(n,l,nu,norm,upper);        //construct upper part
                
                if (n==0)                           //lower part is zero;
                    for(int i=0;i<N;i++)
                        lower.push_back(0.0);
                else
                    transform(n2,l2,nu,norm*alpha,lower);    //construct lower part
            }
            else if(kappa>0){
                E=pow(((2.0/pow(b,2))*(2*n+2*kappa+1)+pow(M,2)),0.5);
                n2=n;
                l=kappa;
                l2=l-1;
                alpha=pow((4*n+4*l+2.0),0.5)/(b*(E+M));
                norm=pow((1.0/(1+pow(alpha,2))),0.5);
                transform(n,l,nu,norm,upper);
                transform(n2,l2,nu,norm*alpha,lower);
            }
        }
        
        else {              //negative energy states
            //cout<<"small"<<endl;
            if(kappa<0){
                E=pow((4*n/pow(b,2)+pow(M,2)),0.5);
                //cout<<E<<endl;
                n2=n-1;
                l=-kappa-1;
                l2=l+1;
                alpha=2.0*pow(n,0.5)/(b*(E+M));
                norm=pow((1.0/(1+pow(alpha,2))),0.5);
                transform(n,l,nu,norm*alpha,upper);
                transform(n2,l2,nu,norm,lower);
            }
            else if(kappa>0){
                E=pow(((2.0/pow(b,2))*(2*n+2*kappa+1)+pow(M,2)),0.5);
                //cout<<E<<endl;
                n2=n;
                l=kappa;
                l2=l-1;
                alpha=-pow((4*n+4*l+2.0),0.5)/(b*(E+M));
                norm=pow((1.0/(1+pow(alpha,2))),0.5);
                transform(n,l,nu,norm*alpha,upper);
                transform(n2,l2,nu,norm,lower);
            }
            
        }
            
        
        E=E*p_sign;
        
    }
    
    
    
};
