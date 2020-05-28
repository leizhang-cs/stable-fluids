#ifndef __FLUID_H__
#define __FLUID_H__

#include "vec.h"
#include <vector>


template<class T, int n>
class Fluid{
    using vec = Vec<T,n>;
    using ivec = Vec<int, n>;
    using pixel = unsigned int;
public:
    Fluid(){}
    Fluid(T visc_input, T kS_input, T aS_input, T dt_input, T density_input,
        vec& L_input, ivec& N_input):
        visc(visc_input), kS(kS_input), aS(aS_input), dt(dt_input), density(density_input),
        L(L_input), N(N_input), size(N_input),
        U0(N_input[0]*N_input[1]), U1(N_input[0]*N_input[1]),
        S0(N_input[0]*N_input[1]), S1(N_input[0]*N_input[1]),
        P(N_input[0]*N_input[1]), div(N_input[0]*N_input[1])
    {
        image_color = new pixel [size[0]*size[1]]();
        for(int i=0; i<n; i++){
            D[i] = L[i]/size[i];
        }
        for(int i=0; i<N[0]; i++){
            for(int j=0; j<N[1]; j++){
                U0[Idx(i,j)][0] = 0;
                U0[Idx(i,j)][1] = 0;
            }
        }
    }

    ~Fluid()
    {
        delete [] image_color;
    }

    pixel* image_color;
    void simulate(vec& F, T Source, vec& X);
    // image as source
    void AddSource(const char* filename);
    void display();
private:
    T visc;
    T kS; // diffusion constant
    T aS; // dissipation constant
    T dt;
    T density;
    // grids row*col: N[0]*N[1]
    vec L, D, O; // D[i] = L[i]/N[i]
    ivec N;
    ivec size;
    std::vector<vec> U0, U1;
    std::vector<T> S0, S1;
    std::vector<T> P; // pressure
    std::vector<T> div; // divergence of velocity

    bool pbc = true; // periodic boundry

    void Vstep(vec& F, T Source, vec& X);
    void AddForce(vec Force, T Source, vec X);
    void Advect();
    void Diffuse();
    void Project();
    void TraceParticle(vec& X1, vec& X0);
    vec Interpolate(std::vector<vec>& U, vec& X);
    // util
    void boundry_condition(std::vector<T>& var);
    void boundry_condition(std::vector<vec>& var);
    int Idx(int i, int j);
    int XtoIdx(vec& X);
};

#endif
