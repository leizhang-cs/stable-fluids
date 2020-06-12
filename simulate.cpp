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
    
    Vec<Type,Ndim> L(256,256);
    Vec<int,Ndim> N(256,256); // D[i] = L[i]/N[i]
    Type visc = 9*1e-4;
    Type kS = 0.1; // diffusion constant
    Type aS = 0.01; // dissipation constant
    Type dt = 1;
    Type density = 128.0;
    //Fluid<Type,Ndim> fluid(visc, kS, aS, dt, density, L, N);
    Vec<Type,Ndim> ForceX(10, 0), ForceY(0, 10), Force(ForceX+ForceY), Force0(0, 0);
    Vec<Type,Ndim> X(0.5, 0.5);
    Type Source = 128, Source0 = 0;
    int simulating;
    std::cout<<"input demo number, 0 for quit: "; std::cin>>simulating;

    using namespace std::this_thread;  // sleep_for, sleep_until
    using namespace std::chrono_literals;

    while(simulating>0){
        Fluid<Type,Ndim> fluid(visc, kS, aS, dt, density, L, N);
        Demo demo = static_cast<Demo>(simulating-1);
        simulating = 1;
        for(int i=0; i<10000; i++){
            /* handle display and user interaction */
            /* get forces F and sources Ssource from the UI */
            switch (demo)
            {
            case Demo::horizontal:
                fluid.simulate(ForceX, Source, X, demo);
                break;
            case Demo::inclined:
                fluid.simulate(Force, Source, X, demo);
                break;
            case Demo::opposite:
                {if(i==0) fluid.AddSource("bunny256.png");
                else fluid.simulate(ForceX, Source0, X, demo);}
                break;
            case Demo::obstacle:
                fluid.simulate(ForceX, Source, X, demo);
                break;
            default:
                fluid.simulate(Force0, Source, X, demo);
                break;
            }
            /*
                if(i==0) fluid.AddSource("bunny128.png");
                else if(true) fluid.simulate(Force, Source0, X);
                else fluid.simulate(Force0, Source0, X);
            */
            if(i%10==0){
                std::cout<<"input 1 for continue; 0 for break: ";
                std::cin>>simulating;
            }
            if(!simulating) break;
            //sleep_for(250ms);
        }
        std::cout<<"input demo number, 0 for quit: "; std::cin>>simulating;
    }

}