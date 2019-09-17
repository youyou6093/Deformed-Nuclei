# Explanation of my functions


## preprocess.h
```c

/* prepare all the hash tables, parameters */
void preprocessing()

/*initial guess of potentials*/
void preprocessing_2(vector<vector<double>> & Phi,vector<vector<double>> &W,vector<vector<double>> &B,vector<vector<double>> & A)

/*generate the scalar and vector potential given other potentials*/
void generate_potential(vector<vector<double>> &Phi,vector<vector<double>> &W, vector<vector<double>> &B,vector<vector<double>> &A,vector<double> &Potential,vector<vector<double>> &scalar_p,vector<vector<double>> &vector_p,int L,int particle_type)
```

## dirac_oscillator.h

```c
dirac_oscillator(int sign,int n,int kappa,double b,double rmin,double rmax,int N)
int E, vector<double> upper, vector<double>lower
double M = 939.0/197.326
```

## states.cpp
```c
State(int a,int b,double c,int d):n(a),k(b),m(c),sign(d)
int n,int k,int sign,double m
string k
```

## generate_state.h
```c
vector<State> generate_statesk(int nodes,int k,double m) 
//I choose kappa max to be 7+max_L
vector<State> generate_statesm(double m,int max_L,int nodes=10)
```


## integrate.h
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


## generate_matrix.h

```c++
/* give m, k1, k2, L,get an key as the input of hash table */
string generate_key(double m,int k1,int k2,int L);

struct params{
int n1,n2,k1,k2,s1,s2;
vector<double>g1,g2,f1,f2;
double m,diag_energy;
vector<vector<double>> scalar_p,vector_p;
};   // a struct that stores the params about the position of the matrix

/* compute the matrix element for 1 channel, 1 position */
double calculate_matrix_element(struct params & my_params,int L)

/* compute the matrix element for 1 position all channels */
void generate_matrix(vector<vector<double>> &M,vector<vector<double>> &scalar_potential,vector<vector<double>> & vector_potential,vector<State> &states,double m,int i,int j)

/* compute the whole matrix */
vector<vector<double>> generate_full_matrix(vector<vector<double>> &scalar_potential,vector<vector<double>> & vector_potential,vector<State> &states,double m)

```




## symm.h

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

bool compare_eig2(eig2 a,eig2 b);

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


## Green-method.h:
```c
/* given density, mass, and pos, solve the klein gordon equation at that point*/
double klein(double mass,vector<double> density,int index);
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

## utility.h:

```c++
/* compute the energy per particle */
double compute_energy(vector<eig2> &occp,vector<eig2> &occn,vector<vector<double>> &Phi,
vector<vector<double>> &W, vector<vector<double>> &B, vector<vector<double>> &A,
vector<vector<double>> &dens, vector<vector<double>> &denv, vector<vector<double>> &den3,
vector<vector<double>> &denp)

/*compute the start point of m */
double magic(int n)

/*flat a 2-d matrix in 1-d , used for matrix diag*/
vector<double> flat_matrix(vector<vector<double>> &M)

/*sort all eigenvalues, git occupied states */
void get_solution(vector<eig2> &occp_raw,vector<eig2> &occn_raw,vector<eig2> &occp,vector<eig2> &occn)

/* for specific m, return the solution eig2,contains m
 basically combine the solution with m*/
vector<eig2> get_temp_solution(vector<eig> &results,double m)

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


## partial result:




# Next step

```
| Model | Nucleus  |  C++  | Python |
|-------|----------|-------|--------|
|Garnet |PB208     |7.892  | 7.890  |
|Garnet |Ca48      |8.615  | 8.614  |
|Garnet |Ca40      |8.530  | 8.529  |
|Garnet |Sn132     |8.362  | 8.359  |
|Gold   |PB208     |7.889  | 7.888  |
|Gold   |Ca48      |8.586  | 8.585  |
|Gold   |Ca40      |8.539  | 8.538  |
|Gold   |Sn132     |8.340  | 8.338  |
|NL3    |PB208     |7.877  | 7.876  |
|NL3    |Ca48      |8.640  | 8.640  |
|NL3    |Ca40      |8.541  | 8.540  |
|NL3    |Sn132     |8.363  | 8.361  |
|Linear |Ca48      |6.712  | 6.712  |
```


## Oct 17
My function can now read the parameters from the a txt file.

## Oct 18
First 50 iterations run the self-consisitent calculations. Then at iteration 51 I add deformation potential

## Oct 19
Try linear parameters

## Oct 24
complete self consistency on Linear conditions, still testing

## Nov 14
Result is not good.

## Nov 22
The calculations of the matrix element is verifyed

## Jun 1st
1. Added parallel

2. Fixed the KG equation and the Possion Equation

3. Current problems, not converge, the shape remains the same




















