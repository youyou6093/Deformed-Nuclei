//tested

#include <iostream>
#include <cmath>
#include <vector>


//THis guy stores coefficients for a Laguerre Polynomial in order to evaluate
//it more efficiently than going recursive.
//As such it will spend a bunch of time in the beginning evaluating all
//coefficients and then do simple lookups (O(1) in time)

//Notes for when I'm old:
//all x^0 coeffs are positive so we can get signs from this
//I'm debating whether it's faster to store actual values or Log of values.
//it's probably the first down the road.
//however they are calculated as exp(log(cc))

//tested with n = 10, a = 31.2 x = 0.22 vs Mathematica function LaguerreL
//result is 1.10592e+09
class LaguerrePolynomialC
{
public:
    int n;//the n value of the polynomial.
    double a;//the alpha value of the polynomial;
    std::vector<double> Coeffs;//stores all coefficients. Coeffs[k] * x^k
    LaguerrePolynomialC():n(-1),a(0){;}//void constructor just in case
    LaguerrePolynomialC(int n_,double a_=0):n(n_),a(a_)
    {
        if (n >= 0)
        {
            Coeffs.resize(n+1,0);
            for (int i = 0; i < Coeffs.size(); ++i)
            {
                Coeffs[i] = std::lgamma(n+a +1.) - std::lgamma(n-i+1.) - std::lgamma(a+i+1.)-std::lgamma(i+1);
                Coeffs[i] = std::exp(Coeffs[i]);
            }
        }
    }
    double Evaluate(double x)
    {
        if (n < 0) return  0;
        double sum = 0, xeval = 1.;
        int phase = 1;
        for (int i = 0; i < Coeffs.size(); ++i)
        {
            sum += phase * Coeffs[i] * xeval;
            xeval *= x;
            phase = -phase;
        }
        return sum;
    }
    //Operator = should work by default as it should copy all members.
};

//Ok this is a bit convoluted.
//when writing this, doing hartree-fock calculations is in the back of my mind
//so each level will NOT contain a copy of the Laguerre polynomial class.
//instead it will be taken as input and all values that are necessary will be stored.
//this might not make much sense the way it written right now because you are not
//in my head so hopefully it will be improved.
//Note to future self h includes end points (counter intuitive given HF calculation)
//Phase is positive near the origin.
class HOLevelWFs
{
public:
    int n;//nodal number (0 for no nodes and so on)
    int l;//angular momentum

    double nu;//HO wave functions go as exp(-nu * r^2)
    double rmin,rmax;
    int N;//number of discretization points.
    double NormalizationFactor;//so that I don't have to calculate it every single time;
    
    std::vector<double> R;
    std::vector<double> dR;//derivative of R (analyticaly calculated.

    HOLevelWFs():n(0),l(0),nu(0.),rmin(0.),rmax(0.),NormalizationFactor(0.){;}
    HOLevelWFs(int n_,int l_,double nu_,double rmin_,double rmax_,int N_)
    :n(n_),l(l_),nu(nu_),rmin(rmin_),rmax(rmax_),N(N_)
    {
        LaguerrePolynomialC Lagpol(n,l+0.5);
        LaguerrePolynomialC Lagpolder(n-1,l+1.5);

        
        NormalizationFactor = 0.5*(std::log(2/M_PI) + 3*std::log(nu)) + (n+2*l+3)*std::log(2);
        NormalizationFactor += std::lgamma(n+1) + l*std::log(nu) - ldifac(2*n+2*l+1);
        NormalizationFactor = std::exp(0.5*NormalizationFactor);
        
        double r = rmin;
        double h = (rmax-rmin)/(N-1.);
//        double fourpi = 4 * M_PI;
        R.resize(N,0);
        dR.resize(N,0);
        
        double L0,L1,nursq;//to store value of Laguerre functions

        
        //Yes, I know this can be changed to something faster but honestly, at this
        //point I'm done worrying about this. It's only going to be calculated once
        //for every single particle level that we're interested in so the time
        //saved is going to be like what? 2 minuites in a 30min calculation?
        //you'd think that, but for N=1024 and 100 evaluations it takes less than 0.1 sec.
        
        for (int i = 0; i < N; ++i)
        {
            nursq = nu*r*r;

            L0 = Lagpol.Evaluate(2*nursq);
            L1 = Lagpolder.Evaluate(2*nursq);
            R[i] = NormalizationFactor * std::pow(r,l) * std::exp(-nursq) * L0;
            
            if (l == 0)
                dR[i] = -NormalizationFactor * std::exp(-nursq) * nu * r *(4*L1 +2 *L0);
            if (l>0)
                dR[i] = -NormalizationFactor * std::exp(-nursq) *std::pow(r,l-1)*(4*nursq*L1 -(l-2*nursq) *L0);

            //This is like horrible. It's like I'm writing in FORTRAN.

            
            r+= h;
            
        }
    }
    //This is here for accessibility
    //not meant to be used really.
    double operator()(double r)
    {
        if (r < rmin || r> rmax)
            std::cerr << "Out of bounds operator call. " << std::endl;
        
        double x = rmin;
        int i = 0;
        double h = (rmax-rmin)/(N-1.);

        while (x < r)
        {
            x+=h;
            ++i;
        }
        --i;//reduce to go where we want. i is now before r and i+1 after it.
        double x1 = rmin + i*h,x2 = x1 + h;
        double y1 = R[i], y2 = R[i+1];
        double slope = (y2-y1)/(x2-x1);
        double offset = y1 - slope * x1;
        
        return slope * r + offset;
    }
private:
    //Guess who crashes after a number too small for our purposes?
    long long difac(long long n)
    {
        if (n==1 || n==0) return 1;
        return n * difac(n-2);
    }
//public:
    double ldifac(int n)
    {
        if (n == 1 || n==0) return 0.;
        return std::log(n) + ldifac(n-2);
    }
    
    
};

//
int test(){
    HOLevelWFs ho(12,8, 1./3., 0., 20., 2048);
    double h = 20./double(2047);
    for (int i = 0; i < ho.N; ++i)
    {
        std::cout << i * h << " " << ho.R[i] << " " << ho.dR[i] << std::endl;
    }
    return 0;
}

