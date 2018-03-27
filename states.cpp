//tested

#ifndef States_h
#define States_h

#include <iostream>
#include <string>
using namespace std;

class State{
public:
    int n,k,sign;
    int m;
    string key;
    State(int a,int b,int c,int d):n(a),k(b),m(c),sign(d){
        key=to_string(n)+"."+to_string(k)+"."+to_string(sign);
    }
    
    State(){
        
    }
    
    friend ostream& operator<<(ostream &os, State x){
        os<<x.n<<' '<<x.k<<' '<<x.m<<' '<<x.sign;
        
        return os;
    }
    
    
    
};


class State_simp{
public:
    int n,k,sign;
    State_simp(int a,int b,int c):n(a),k(b),sign(c){}
};



#endif
