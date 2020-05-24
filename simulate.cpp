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

    //std::cout<<"beforeAlloc"<<std::endl;
    // 32,32. 32, 32. dt=1, density=1, F=(1,0.5)
    Vec<Type,Ndim> L(0.32,0.32);
    Vec<int,Ndim> N(128,128); // D[i] = L[i]/N[i]
    Type visc = 9*1e-4;
    Type kS = 0.01; // diffusion constant
    Type aS = 0.1; // dissipation constant
    Type dt = 0.01;
    Type density = 1000.0;
    //std::cout<<"beforeFluid"<<L[0]<<" "<<L[1]<<" "<<N[0]<<std::endl;
    Fluid<Type,Ndim> fluid(visc, kS, aS, dt, density, L, N);
    Vec<Type,Ndim> Force(-100, 100), F0(0, 0);
    Vec<Type,Ndim> X(0.5, 0.5);
    Type Source = 1e4, Source0 = 0;
    int simulating = 1;
    
    std::cout<<"beforeLoop, L:"<<L[0]<<", N:"<<N[0]<<std::endl;
{
    using namespace std::this_thread;     // sleep_for, sleep_until
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
        else if(true) fluid.simulate(Force, Source, X);
        else fluid.simulate(F0, Source0, X);
        
        fluid.display();

        sleep_for(250ms);
    }
}

}