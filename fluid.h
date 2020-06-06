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
        L(L_input), N(N_input),
        U0(N_input[0]*N_input[1]), U1(N_input[0]*N_input[1]),
        S0(N_input[0]*N_input[1]), S1(N_input[0]*N_input[1]),
        P(N_input[0]*N_input[1]), div(N_input[0]*N_input[1])
    {
        image_color = new pixel [N[0]*N[1]]();
        for(int i=0; i<n; i++){
            D[i] = L[i]/N[i];
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
    vec L; // grids row*col: N[0]*N[1], TODO: L correspond to X
    vec D; // D[i] = L[i]/N[i]
    ivec N;
    std::vector<vec> U0, U1;
    std::vector<T> S0, S1;
    std::vector<T> P; // pressure
    std::vector<T> div; // divergence of velocity

    bool pbc = true; // periodic boundry
    bool obstacle = false; // obstacle area
    bool debug_flag = false;

    void Vstep(vec& F, T Source, vec& X);
    void AddForce(vec Force, T Source, vec X);
    void Advect();
    void Diffuse();
    void Project();

    void TraceParticle(vec& X1, vec& X0);
    vec Interpolate(std::vector<vec>& U, vec& X);
    T Interpolate(std::vector<T>& S, vec& X);
    void boundary_condition(std::vector<vec>& var); // TODO: aperiodic boundary
    void boundary_condition(std::vector<T>& var); 
    // util
    int Idx(int i, int j);
    int XtoIdx(const vec& X);
    T IdxtoX(const int i, const int dim); // i,j to y,x
    void Check_Symmetry(std::vector<vec>& U, std::vector<T>& S, const std::string& s);
};

#endif
