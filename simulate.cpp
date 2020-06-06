#include "fluid.h"
#include "vec.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <thread>


int main(){
    using Type = double;
    const int Ndim = 2;
    
    Vec<Type,Ndim> L(128,128);
    Vec<int,Ndim> N(64,128); // D[i] = L[i]/N[i]
    Type visc = 9*1e-4;
    Type kS = 0.1; // diffusion constant
    Type aS = 0.01; // dissipation constant
    Type dt = 1;
    Type density = 1000.0;
    Fluid<Type,Ndim> fluid(visc, kS, aS, dt, density, L, N);
    Vec<Type,Ndim> Force(-0.2, 0), Force0(0, 0);
    Vec<Type,Ndim> X(0.5, 0.5);
    Type Source = 1024, Source0 = 0;
    int simulating = 1;
    
    std::cout<<"beforeLoop, L:"<<L[0]<<", N:"<<N[0]<<std::endl;
{
    using namespace std::this_thread;  // sleep_for, sleep_until
    using namespace std::chrono_literals;
    for(int i=0; i<10000; i++){
        /* handle display and user interaction */
        /* get forces F and sources Ssource from the UI */
        if(i%10==0){
            std::cout<<"input 1 for continue; 0 for break:";
            std::cin>>simulating;
        }
        if(!simulating) break;
        
        if(i==0) fluid.AddSource("bunny128.png");
        else if(true) fluid.simulate(Force, Source0, X);
        else fluid.simulate(Force0, Source0, X);
        
        fluid.display();

        sleep_for(250ms);
    }
}

}