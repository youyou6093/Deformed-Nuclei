//tested

#ifndef my_spline_h
#define my_spline_h

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <iostream>
using namespace std;

class my_spline{
    
private:
    //spline
    int size;   //size of y array and x array
    double *y,*x;
    double tolerance;
    //integral
    double result, error;
    struct sp{
        gsl_interp_accel *acc;
        gsl_spline *spline;
    }sp_object;
    
    
    
    
    
    void start_spline(){            //set up the spline
        sp_object.acc = gsl_interp_accel_alloc ();
        sp_object.spline = gsl_spline_alloc (gsl_interp_cspline,size);
        gsl_spline_init (sp_object.spline, x, y, size);
        //now the initialization of spline is finished;
    }
    
    double p_integral(double range1,double range2){
    	gsl_set_error_handler_off();   // this one is pretty important
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
        // this 1000 might need to change;
        gsl_function F;
        F.params=&sp_object;
        F.function=eval2;
        gsl_integration_qags (&F, range1, range2, 0, tolerance, 10000, w, &result, &error);
        //the acuraccy
        //this 1000 might need to be changed;
        gsl_integration_workspace_free(w);
        return result;             //get the result of an integration between range1 and 2
    }
    
    
public:
    //constructor from vectors
    my_spline(vector<double> &y_vector,vector<double> &x_vector,double tol){
        size=int(x_vector.size());
        y=&y_vector[0];
        x=&x_vector[0];
        start_spline();
        tolerance=tol;
    }
    
    //constructor from arrays
    my_spline(double *y_in,double *x_in,int length,double tol){
        y=y_in;
        x=x_in;
        size=length;
        start_spline();
        tolerance=tol;
    }
    
    double eval(double value){
        return gsl_spline_eval (sp_object.spline, value, sp_object.acc);
    }    //this is not I am using for spline;
    
    
    //this is the function I put into the integration
    static double eval2(double value,void *params){
        sp * v = (sp*) params;        //cast params into sp*
        return gsl_spline_eval (v->spline, value, v->acc);
    }
    
    
    ~my_spline(){
        gsl_spline_free (sp_object.spline);
        gsl_interp_accel_free (sp_object.acc);
        
        //free memory
    }
    
    double integral(double range1,double range2){
        return p_integral(range1, range2);
    }
    
    double integral(){
        return p_integral(x[0], x[size-1]);
    }
    
    
    
    
};

#endif


