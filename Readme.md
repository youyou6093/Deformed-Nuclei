# Explanation of my functions


## dirac_oscillator

```c
dirac_oscillator(int sign,int n,int kappa,double b,double rmin,double rmax,int N)
int E, vector<double> upper, vector<double>lower
double M = 939.0/197.326
```

## State
```c
State(int a,int b,double c,int d):n(a),k(b),m(c),sign(d)
int n,int k,int sign,double m
string k
```

## Generate_state.h
```c
vector<State> generate_statesk(int nodes,int k,double m) 
//I choose kappa max to be 7+max_L
vector<State> generate_statesm(double m,int max_L,int nodes=10)
```


## my_spline(integrate.h)
```c
my_spline(vector<double> &y_vector,vector<double> &x_vector,double tol)
my_spline(double *y_in,double *x_in,int length,double tol)
double integral(double range1,double range2)
double integral() //integrate over the range of x
```


## simps.h
```c
double simps(vector<double> &y,vector<double> &x)  //size of the array has to be odd
```

## matrix_diag(symm.h)

```c
struct eig{              //a structure that stores the eigenvectors
    double eigen_values;
    vector<double> eigen_vectors;
};
vector<eig> results

matrix_diag(double *data,double n)    
matrix_diag(vector<double> &data,double n)
matrix_diag(vector<vector<double>> &M,int size)
get_results()  //doesn’t return anything
```

## Solution.h

```c
struct eig2{
    eig solution;
    double m;
};

struct coef_pair{             //takes coefficient and state
    double coefs;
    State state;
};

struct solution_wave_function{                // a struct that stores wave function for one kappa
    int kappa;
    vector<double>upper;
    vector<double>lower;
};
```

### `class Solution` (Get the wave functions from eigenvectors)

```c++
class Solution{
public:
    vector<coef_pair> my_pair;
    double m;     //the quantum number of this solution
    double energy;
    vector<int> kappas;
    State primary_state;
    vector<solution_wave_function> wavefunctions;

Solution(eig2 result)
bool add_kappa(int k)         //decide whether to add new kappa into collection
void get_primary_state()      //find the state with the biggest coefficients
solution_wave_function get_wave_function(int kappa) //get wavefunction with respect to one kappa
void get_all_wave_function()        //get wavefunctions with all possible kappa
}
```




## Density.h

```c
vector<Solution> get_solutions_object(vector<eig2> occ)    //get all the wavefunctions based on the eigenvectos
vector<int> get_possible_L(int k1,int k2)                  
//give kappa1 and kappa2 , get all possible L, didn't consider parity
struct den{                   //density for one channel
    vector<double> s;
    vector<double> v;
    //int channel;
};
void generate_density(vector<eig2> &occn, vector<eig2> &occp, vector<vector<double>> &dens, vector<vector<double>> &denv,vector<vector<double>> &denp,vector<vector<double>> &den3)
//given occs, compute the densities

```





### `class Density`

```c
class Density{
public:
    vector<Solution> solution;
    vector<den> density;
    Density(vector<eig2> & occ);
    //add the new desnity calculated form 1 single occupied state
    void append(den temp,int channel);
    void compute_one_(int a,int b,int num)； 
    //compute a partial density according to a specific soluiton and two specific kappas in that solution.
    void compute(int num)； //compute the density for 1 occ state.
    void Compute_all();     //compute the total densities
```



	

## effecitve_density.h:
```c
extern"C" {
    double tj_(double* a1, double* a2, double* a3,double* b1, double* b2, double* b3);
    double sj_(double* j1, double* j2, double* j12,double* j3,double* j,double* j23);
    double coef9_(double* a1,double* a2,double* a3,
                  double* a4,double* a5,double* a6,
                  double* a7,double* a8,double* a9);
}         //declare that there is something in the .o library

double pcg_(double j1,double j2,double j3) //cg coefficients
vector<int> possible_L(int l1,int l2) //get all possible l based on l1,l2
vector<vector<double>> multiplication(vector<vector<double>> &vec1,vector<vector<double>> &vec2)
//multiply two lengendre series
void get_effective_density(vector<vector<double>> &EFF_Phi,vector<vector<double>>  &EFF_B,vector<vector<double>> &EFF_A,vector<vector<double>> &EFF_W, vector<vector<double>> &Phi,vector<vector<double>> &W,vector<vector<double>> &B,vector<vector<double>> &A, vector<vector<double>> &dens,vector<vector<double>> &denv,vector<vector<double>> & den3,vector<vector<double>> denp)
//Get effective density based on the non-linear part
void get_potential(vector<vector<double>> &EFF_Phi,vector<vector<double>>  &EFF_B,vector<vector<double>> &EFF_A,vector<vector<double>> &EFF_W, vector<vector<double>> &Phi,vector<vector<double>> &W,vector<vector<double>> &B,vector<vector<double>> &A, vector<vector<double>> &dens,vector<vector<double>> &denv,vector<vector<double>> & den3,vector<vector<double>> denp, int L)
//update the potential based on the density
void update_potential(...)     // a small iteration, keep compute the densities and get potential

```




## main.cpp:


```c
unordered_map<string, vector<vector<double>> > States_map;   //map from states to upper,lower pair
unordered_map<string, double> Energy_map;                    //map from states to energy
unordered_map<string, double > Angular_map;                  //map from m,k1,k2,L to angular term

void preprocessing(); 										//construct all the different maps
void preprocessing_2(vector<vector<double>> & Phi,vector<vector<double>> &W,vector<vector<double>> &B,vector<vector<double>> & A);
// initialize all the potentials

void generate_potential(vector<vector<double>> &Phi,vector<vector<double>> &W, vector<vector<double>> &B,vector<vector<double>> &A,vector<double> &Potential,vector<vector<double>> &scalar_p,vector<vector<double>> &vector_p,int L,int particle_type) //L represent the potential channel
//construct the scalar and vector potentials
```


# Next step

#### Let the function read the parameters
#### add threads
#### find the reason why the program doesn't work on deformed condition
