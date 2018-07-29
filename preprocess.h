//
//  preprocess.h
//  
//
//  Created by Junjie Yang on 10/25/17.
#ifndef preprocess_h
#define preprocess_h

/*construct three unordered map: States, energy and Angular and get global data */
void preprocessing(){
    /* get the input parameters */
    ifstream parameter_file;
    parameter_file.open(parameters);
    parameter_file >> ms;
    parameter_file >> mv;
    parameter_file >> mp;
    parameter_file >> mg;
    parameter_file >> gs;
    parameter_file >> gv;
    parameter_file >> gp;
    parameter_file >> gg;
    parameter_file >> lambdas;
    parameter_file >> lambdav;
    parameter_file >> lambda;
    parameter_file >> ks;
    parameter_file >> ka;
    /*--------------------------------------------- */
    cout << ms << ' ' << mv << ' ' << mp << ' ' << mg << endl;
    cout << gs << ' ' << gv << ' ' << gp << ' ' << gg << endl;
    cout<< lambdas << ' ' << lambdav << ' ' << lambda << ' ' << ks << ' ' << ka << endl;
    ms /= hbarc;
    mv /= hbarc;
    mp /= hbarc;
    mg /= hbarc;
    ka /= hbarc;
    
    
    vector<vector<double>> wave_function;  //store the upper_part and lower part
    vector<double> energy;                 //stor all the energys
    for(int i=0;i<N;i++)
        fx.push_back(h*i+rmin);   //construct the x array for all the basis
    
    /*constructed the States hash table*/
    
    /*this is the largest basis, gives me everything in the basis,
      MAX_KAPPA is defined in the main.cpp*/
    vector<State>states_1 = generate_statesm(1,MAX_KAPPA);     //right now it is 7
    for(int i=0;i<states_1.size();i++){
        dirac_oscillator ho(states_1[i].sign,states_1[i].n,states_1[i].k,b,rmin,rmax,N);
        wave_function.clear();   //empty the array
        wave_function.push_back(ho.upper);
        wave_function.push_back(ho.lower);
        States_map.insert(make_pair(states_1[i].key, wave_function));
        Energy_map.insert(make_pair(states_1[i].key, ho.E));
    }
}



/*initialize the potentials to wood-saxon form*/
void preprocessing_2(vector<vector<double>> & Phi,vector<vector<double>> &W,vector<vector<double>> &B,vector<vector<double>> & A){
    double value;
    
    /* wood saxon part */
    for(int i=0;i<N;i++){
        value=(exp(fx[i]-4.3)/0.7+1);
        Phi[0][i] = start[0] / value;
        W[0][i] = start[1] / value;
        B[0][i] = start[2] / value;
        A[0][i] = start[3] / value;
    }
}



/*generate the scalar and vector potentials based on PHI W B A
 L is the channel of the potential we add*/
void generate_potential(vector<vector<double>> &Phi,vector<vector<double>> &W, vector<vector<double>> &B,vector<vector<double>> &A,vector<double> &Potential,vector<vector<double>> &scalar_p,vector<vector<double>> &vector_p,int L,int particle_type){
    

    if (particle_type==0){                                //neutron
        for(int i=0; i<max_L; i++){                    //for all the channel
            for(int j=0; j<N; j++){
            	scalar_p[i][j] = -Phi[i][j];
            	vector_p[i][j] = W[i][j]-0.5*B[i][j];
            }
        }
    }
    else{                                                //proton
        for(int i=0;i<max_L;i++){
            for(int j=0;j<N;j++){
            	scalar_p[i][j] = -Phi[i][j];
            	vector_p[i][j] = W[i][j]+0.5*B[i][j]+A[i][j];
            }
        }
    }
    
    /* USE vector potential*/
    /* put the potential we add to the scalar and vector part, depend on what I want the potential to be*/
    if (L < max_L){
        if(particle_type == 0){
            for(int i=0;i<N;i++){
               vector_p[L][i] -= Potential[i];
            }
        }
        else{
            for( int i =0; i < N; i++){
                vector_p[L][i] += Potential[i];
            }
        }
    }
}



#endif /* preprocess_h */
