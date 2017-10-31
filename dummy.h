//
//  dummy.h
//  
//
//  Created by Junjie Yang on 10/25/17.
//

#ifndef dummy_h
#define dummy_h


/* main */
/* output potential */
//ofstream myfile;
//myfile.open("potential.txt");
//for(int i=0;i<N;i++)
//    myfile<<fx[i]<<' '<<Phi[0][i]<<' '<<W[0][i]<<' '<<B[0][i]<<' '<<A[0][i]<<endl;


/* output density */
//        ofstream myfile;
//        myfile.open("density.txt");
//        for(int i=0;i<N;i++){
//            myfile<<fx[i]<<' '<<dens[0][i]<<' '<<denv[0][i]<<' '<<den3[0][i]<<' '<<denp[0][i]<<endl;
//        }
//        myfile.close();

/* output the single particle energy */
//cout << "proton_energy" << endl;
//for(int i=0;i<occp.size();i++) cout<<Final_occp[i].energy<<' '<<Final_occp[i].m<<endl;
//cout << "neutron_energy" << endl;
//for(int i=0;i<occn.size();i++) cout<<Final_occn[i].energy<<' '<<Final_occn[i].m<<endl;
/* ----------------------------------- */


//        /*test*/
//        Solution test = Final_occp[0];
//        test.get_primary_state();
//        cout << test.energy << ' ' << test.m << ' ' << test.primary_state << endl;
//        test.get_all_wavefunction();
//        for( int i =0; i < test.my_pair.size(); i++)
//            cout << test.my_pair[i].coefs << ' ' << test.my_pair[i].state << endl;
//        for( int i = 0; i < test.wavefunctions.size(); i++){
//            ofstream infile;
//            infile.open(to_string(test.wavefunctions[i].kappa) + "wavefunctionc.txt");
//            for( int x = 0 ; x < N; x++){
//                infile << fx[x] << ' ' << test.wavefunctions[i].upper[x] << ' ' << test.wavefunctions[i].lower[x] << endl;
//            }
//            infile.close();
//        }
////        for(int i = 0; i < N; i++){
////            cout << scalar_n[0][i] << ' ' << vector_n[0][i] << ' ' << scalar_p[0][i] << ' ' << vector_p[0][i] << endl;
////        }
//        break;
//





#endif /* dummy_h */
