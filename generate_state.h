//tested

#ifndef generate_state_h
#define generate_state_h



#include<vector>
#include "states.cpp"
using namespace std;

vector<State> generate_statesk(int nodes,int k,double m){
    vector<State> negative,positive,my_state;
    if (k>0)
        negative.push_back(State(0,k,m,-1));
    else
        negative.push_back(State(1,k,m,-1));
    positive.push_back(State(0,k,m,1));
    for(int i=0;i<nodes;i++){
        positive.push_back(State(positive[i].n+1,k,m,positive[i].sign));
        negative.push_back(State(negative[i].n+1,k,m,negative[i].sign));
    }
    my_state=positive;
    my_state.insert(my_state.end(),negative.begin(),negative.end());
    //cout<<k<<' '<<my_state.size()<<endl;
    return my_state;

    
}

vector<State> generate_statesm(double m,int max_L,int nodes=10){
    vector<State> my_states,state1;
    int kmin=int(abs(m)+0.5);
    int kmax=7+max_L;
    for(int i=kmin;i<kmax+1;i++){
        state1=generate_statesk(nodes,i,m);   //i is k
        my_states.insert(my_states.end(),state1.begin(),state1.end());
        state1=generate_statesk(nodes,-i,m);  //-i is another k
        my_states.insert(my_states.end(),state1.begin(),state1.end());
    }
    return my_states;
};


#endif
