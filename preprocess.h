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
//    double something;
    parameter_file.open(parameters);
    parameter_file >> ms;
//    cout << something << endl;
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
    vector<State>states_1 = generate_statesm(1,MAX_KAPPA);  
    for(int i=0;i<states_1.size();i++){
        dirac_oscillator ho(states_1[i].sign,states_1[i].n,states_1[i].k,b,rmin,rmax,N);
        wave_function.clear();   //empty the array
        wave_function.push_back(ho.upper);
        wave_function.push_back(ho.lower);
        States_map.insert(make_pair(states_1[i].key, wave_function));
        Energy_map.insert(make_pair(states_1[i].key, ho.E));
    }
    // /*construct the angular hash table, basically read values from the file and store it inside unordered set*/
    // vector<string> keys_1,value_1;
    // string word;
    // ifstream infile;
    // infile.open("test.txt");
    // while(infile>>word){
    //     keys_1.push_back(word);
    //     infile>>word;
    //     value_1.push_back(word);
    // }
    // infile.close();
    // for(int i=0;i<keys_1.size();i++){
    //     Angular_map.insert(make_pair(keys_1[i],stod(value_1[i])));
    // }
}



/*initialize the potentials to wood-saxon form*/
void preprocessing_2(vector<vector<double>> & Phi,vector<vector<double>> &W,vector<vector<double>> &B,vector<vector<double>> & A){
    double value;
    
    /* clear the Potentials */
    Phi.clear();
    W.clear();
    B.clear();
    A.clear();
    vector<double> phi,w,b,a;
    /* wood saxon part */
    for(int i=0;i<N;i++){
        value=(exp(fx[i]-4.3)/0.7+1);
        phi.push_back(start[0]/value);
        w.push_back(start[1]/value);
        b.push_back(start[2]/value);
        a.push_back(start[3]/value);
    }
    Phi.push_back(phi);
    W.push_back(w);
    B.push_back(b);
    A.push_back(a);
    phi.clear();
    w.clear();
    b.clear();
    a.clear();
    /*initialize the 0 channel to wood saxon form*/
    
    /*initialize all other channels to 0*/
    for(int i=0;i<N;i++)
        empty.push_back(0.0);
    for(int i=0;i<max_L-1;i++){         //max_L - 1 because I already assign the 0 channel
        Phi.push_back(empty);
        W.push_back(empty);
        B.push_back(empty);
        A.push_back(empty);
    }
    
    
    /*for test purpose, make an guess at the beginning*/
//    for(int i = 0; i < N; i ++){
//        Phi[1][i] = -200 * Deformation_parameter * fx[i] * fx[i] * exp(-fx[i]);
//    }
}



/*generate the scalar and vector potentials based on PHI W B A
 L is the channel of the potential we add*/
void generate_potential(vector<vector<double>> &Phi,vector<vector<double>> &W, vector<vector<double>> &B,vector<vector<double>> &A,vector<double> &Potential,vector<vector<double>> &scalar_p,vector<vector<double>> &vector_p,int L,int particle_type){
    
    /*clear the scalar potential and vector potential at the beginning, very important*/
    scalar_p.clear();
    vector_p.clear();
    vector<double> temp1,temp2;
    if (particle_type==0){                                //neutron
        for(int i=0;i<Phi.size();i++){                    //for all the channel
            /*temp1 and temp2 is the potential for 1 channel*/
            temp1.clear();                                //clear the container
            temp2.clear();
            for(int j=0;j<Phi[0].size();j++){
                temp1.push_back(-Phi[i][j]);
                temp2.push_back(W[i][j]-0.5*B[i][j]);
            }
            scalar_p.push_back(temp1);
            vector_p.push_back(temp2);
        }
    }
    else{                                                //proton
        for(int i=0;i<Phi.size();i++){
            temp1.clear();
            temp2.clear();
            for(int j=0;j<Phi[0].size();j++){
                temp1.push_back(-Phi[i][j]);
                temp2.push_back(W[i][j]+0.5*B[i][j]+A[i][j]);
            }
            scalar_p.push_back(temp1);
            vector_p.push_back(temp2);
        }
    }
    
    /* USE vector potential*/
    /* put the potential we add to the scalar and vector part, depend on what I want the potential to be*/
    if (L < max_L){
        if(particle_type == 0){
            for(int i=0;i<N;i++){
               vector_p[L][i] -= Potential[i];
                // vector_p[L][i] += 0;
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
