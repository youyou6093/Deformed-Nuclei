/* I think I finished*/


#ifndef generate_state_h
#define generate_state_h



#include<vector>
#include "states.cpp"
using namespace std;
extern int max_k;


/*generate the states with the same k*/
vector<State> generate_statesk(int nodes,int k,int m){
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
    return my_state;

    
}

vector<State> generate_statesm(int m,int max_L,int nodes=10){
    vector<State> my_states,state1;
    int kmin = (abs(m) + 1) / 2;
    int kmax = max_k + max_L;
    for(int i=kmin;i<kmax+1;i++){
        /* For every k I choose, there is positive k and negative k*/
        state1=generate_statesk(nodes,i,m);
        my_states.insert(my_states.end(),state1.begin(),state1.end());
        state1=generate_statesk(nodes,-i,m);
        my_states.insert(my_states.end(),state1.begin(),state1.end());
    }
    return my_states;
};


#endif
