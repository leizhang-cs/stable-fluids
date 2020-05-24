#include "fluid.h"
#include "vec.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>


int main(){
    using Type = double;
    const int Ndim = 2;

    //std::cout<<"beforeAlloc"<<std::endl;
    // 32,32. 32, 32. dt=1, density=1, F=(1,0.5)
    Vec<Type,Ndim> L(320,320), N(32,32); // D[i] = L[i]/N[i]
    Type visc = 5.0f;
    Type kS = 0; // diffusion constant
    Type aS = 0; // dissipation constant
    Type dt = 0.001f;
    Type density = 1.0f;
    //std::cout<<"beforeFluid"<<L[0]<<" "<<L[1]<<" "<<N[0]<<std::endl;
    Fluid<Type,Ndim> fluid(visc, kS, aS, dt, density, L, N);
    Vec<Type,Ndim> Force(0.1, 0.05), F0(0, 0);
    Vec<Type,Ndim> X(0.5, 0.5);
    Type Source = 0;
    int simulating = 1;
    
    std::cout<<"beforeLoop, L:"<<L[0]<<", N:"<<N[0]<<std::endl;
    for(int i=0; simulating && i<10000; i++){
        /* handle display and user interaction */
        /* get forces F and sources Ssource from the UI */
        if(i%10==0){
            std::cout<<"input 1 for continue; 0 for break:";
            std::cin>>simulating;
        }

        fluid.simulate(Force, Source, X);
        //fluid.simulate(F0, Source, X);
        fluid.display();
    }
}