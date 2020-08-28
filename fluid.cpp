#include "fluid.h"
#include "pixel.h"
#include <fftw3.h>
#include <cassert>


using Type = double;
const int Ndim = 2;

template class Fluid<Type,Ndim>;
using vec = Vec<Type,Ndim>;
using pixel = unsigned int;

static const double small_t = 1e-4;


template<class T, int n>
void Fluid<T,n>::simulate(vec& F, T Source, vec& X, Demo d){
    demo = d;
    if(demo==Demo::obstacle) obstacle = true;
    Vstep(F, Source, X);
    std::swap(U0, U1);
    std::swap(S0, S1);
    display();
    obstacle = false;
}


template<class T, int n>
void Fluid<T,n>::Vstep(vec& F, T S, vec& X){
    // U0, U1, visc, F, dt
    AddForce(F, S, X);
    Advect();
    Diffuse();
    Project();
}


template<class T, int n>
void Fluid<T,n>::AddForce(vec F, T S, vec X){
    // U1 = U0 + F*dt
    // initialize
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            U1[Idx(i,j)].make_zero();
            S1[Idx(i,j)] = 0;
        }
    }

    // w1 = f(w0)
    vec du_top = F*dt;
    for(int i=0; i<N[0]/3; i++){
        for(int j=0; j<N[1]; j++){
            U1[Idx(i,j)] = du_top/S0[Idx(i,j)];
        }
    }

    if(demo!=Demo::opposite){
        for(int i=N[0]/3; i<2*N[0]/3; i++){
            for(int j=0; j<N[1]; j++){
                U1[Idx(i,j)] = du_top/S0[Idx(i,j)];
            }
        }
    }
    
    vec du_btm = demo==Demo::opposite? -du_top: du_top;
    for(int i=2*N[0]/3; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            U1[Idx(i,j)] = du_btm/S0[Idx(i,j)];
        }
    }
    
    // add source
    T ds = S;
    if(demo==Demo::inclined){
        for(int i=3*N[0]/4; i<N[0]; i++){
            for(int j=0; j<N[1]/4; j++){
                S1[Idx(i,j)] = ds;
            }
        }
    }
    else{
        for(int i=3*N[0]/8; i<5*N[0]/8; i++){
            for(int j=0; j<N[1]/4; j++){
                S1[Idx(i,j)] = ds;
            }
        }
    }
    
    // obstacle
    if(obstacle){
        for(int i=0; i<N[0]; i++)
            for(int j=0; j<N[1]; j++)
                if((i-N[0]/2)*(i-N[0]/2)+(j-N[1]/2)*(j-N[1]/2)<std::min(N[0],N[1])){
                    U0[Idx(i,j)] = vec(0,0);   
                }
    }
    //Check_Symmetry(U1, S1, "addF");
    std::cout<<"du top: "<<du_top<<", btm: "<<du_btm<<std::endl;
    std::cout<<"delta source: "<<ds<<std::endl;
}

template<class T, int n>
void Fluid<T,n>::Advect(){
    // U1, U0, dt
    // TraceParticle: method of charactristic
    std::vector<vec> Ut(N[0]*N[1]);

    for(int i=0; i<N[0]; i++){
        T y = IdxtoX(i,0); // (0, N[1]) <=> (L[1], -L[1])
        for(int j=0; j<N[1]; j++){
            vec X0;
            vec X1(IdxtoX(j,1), y); // (0, N[0]) <=> (-L[0], L[0])
            
            TraceParticle(X1, X0);

            // w2 = f(w1)
            // Interpolation
            U1[Idx(i,j)] += Interpolate(U0, X0);
            S1[Idx(i,j)] += Interpolate(S0, X0);
        }
    }
    //Check_Symmetry(U1, S1, "adv");
    std::cout<<"U After Adv, top:"<<U1[Idx(N[0]/3-10,N[1]/2+10)]<<",btm:"<<U1[2*Idx(N[0]/3+10,N[1]/2)]<<std::endl;
}


template<class T, int n>
void Fluid<T,n>::Diffuse(){
    // U1, U0, visc, dt 
    // FTCS, BTCS, FFT

    if(!FFT_scheme){
        // FTCS scheme
        // T k = visc*dt;
        // std::vector<vec> Ut(U1);
        // for(int i=0; i<N[0]; i++){
        //     for(int j=0; j<N[1]; j++){
        //         // TODO: optimization
        //         T ux = (Ut[Idx(i+1,j)][0] - 2.0*Ut[Idx(i,j)][0] + Ut[Idx(i-1,j)][0])/(D[0]*D[0]);
        //         T uy = (Ut[Idx(i,j+1)][1] - 2.0*Ut[Idx(i,j)][1] + Ut[Idx(i,j-1)][1])/(D[1]*D[1]);
        //         U1[Idx(i,j)] += (1+k)*vec(ux,uy);
        //     }
        // }

        // BTCS scheme
        T k = visc*dt;
        std::vector<vec> Ut(U1);
        for(int iter=0; iter<20; iter++){
            for(int i=0; i<N[0]; i++){
                for(int j=0; j<N[1]; j++){
                    U1[Idx(i,j)][0] = Ut[Idx(i,j)][0] + 
                        k * (U1[Idx(i+1,j)][0] - 2.0*U1[Idx(i,j)][0] + U1[Idx(i-1,j)][0])/(D[0]*D[0]);
                    U1[Idx(i,j)][1] = Ut[Idx(i,j)][1] + 
                        k * (U1[Idx(i,j+1)][1] - 2.0*U1[Idx(i,j)][1] + U1[Idx(i,j-1)][1])/(D[1]*D[1]);
                }
            }
        }
    }
    else{
        int num = N[0]*N[1];
        fftw_complex* ux = new fftw_complex[num]();
        fftw_complex* uy = new fftw_complex[num]();
        fftw_complex* s = new fftw_complex[num]();
        fftw_complex* uxf = new fftw_complex[num]();
        fftw_complex* uyf = new fftw_complex[num]();
        fftw_complex* sf = new fftw_complex[num]();

        for(int i=0; i<N[0]; i++){
            for(int j=0; j<N[1]; j++){
                ux[Idx(i,j)][0] = U1[Idx(i,j)][0];
                uy[Idx(i,j)][0] = U1[Idx(i,j)][1];
                s[Idx(i,j)][0] = S1[Idx(i,j)];
            }
        }
        //std::cout<<"U before Diff, top:"<<U1[Idx(N[0]/3-10,N[1]/2+10)]<<",btm:"<<U1[Idx(N[0]/3+10,N[1]/2)]<<std::endl;
        //std::cout<<"S before Diff: "<<min_val<<std::endl;
        fftw_plan pf_x = fftw_plan_dft_1d(num, ux, uxf, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan pf_y = fftw_plan_dft_1d(num, uy, uyf, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan pf_s = fftw_plan_dft_1d(num, s, sf, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(pf_x);
        fftw_execute(pf_y);
        fftw_execute(pf_s);

        for(int i=0; i<N[0]; i++){
            //T x = i;
            //T y = IdxtoX(i,0);
            T y = i<=N[0]/2? i: i-N[0];
            for(int j=0; j<N[1]; j++){
                // w(k) - k*w(k)*k
                //T x = IdxtoX(j,1);
                T x = 0.5*j;
                T r = x*x + y*y;
                T f = 1/(1 + visc*dt*r);
                uxf[Idx(i,j)][0] *= f; uxf[Idx(i,j)][1] *= f;
                uyf[Idx(i,j)][0] *= f; uyf[Idx(i,j)][1] *= f;
            }
        }// Diffuse
        
        //Projection by FFT
        // for(int i=0; i<N[0]; i++){
        //     //T x = i;
        //     //T y = IdxtoX(i,0);
        //     T y = i<=N[0]/2? i: i-N[0];
        //     for(int j=0; j<N[1]; j++){
        //         // w(k) - k*w(k)*k
        //         //T x = IdxtoX(j,1);
        //         T x = 0.5*j;
        //         T r = x*x + y*y;
        //         if(r==0.0) continue;
        //         fftw_complex uxt, uyt;
        //         uxt[0] = (1-x*x/r)*uxf[Idx(i,j)][0] - x*y/r*uyf[Idx(i,j)][0];
        //         uxt[1] = (1-x*x/r)*uxf[Idx(i,j)][1] - x*y/r*uyf[Idx(i,j)][1];
        //         uyt[0] = -y*x/r*uxf[Idx(i,j)][0] + (1-y*y/r)*uyf[Idx(i,j)][0];
        //         uyt[1] = -y*x/r*uxf[Idx(i,j)][1] + (1-y*y/r)*uyf[Idx(i,j)][1];
        //         uxf[Idx(i,j)][0] = uxt[0]; 
        //         uxf[Idx(i,j)][1] = uxt[1];
        //         uyf[Idx(i,j)][0] = uyt[0]; 
        //         uyf[Idx(i,j)][1] = uyt[1];
        //     }
        // }
        
        fftw_plan pb_x = fftw_plan_dft_1d(num, uxf, ux, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan pb_y = fftw_plan_dft_1d(num, uyf, uy, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan pb_s = fftw_plan_dft_1d(num, sf, s, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(pb_x);
        fftw_execute(pb_y);
        fftw_execute(pb_s);
        
        for(int i=0; i<N[0]; i++){
            for(int j=0; j<N[1]; j++){
                U1[Idx(i,j)][0] = ux[Idx(i,j)][0]/num;
                U1[Idx(i,j)][1] = uy[Idx(i,j)][0]/num;
                S1[Idx(i,j)] = s[Idx(i,j)][0]/num;
            }
        }    
        std::cout<<"U After Proj, top:"<<U1[Idx(N[0]/3-10,N[1]/2+10)]<<",btm:"<<U1[2*Idx(N[0]/3+10,N[1]/2)]<<std::endl;
        //Check_Symmetry(U1, S1, "diff");

        fftw_destroy_plan(pf_x);
        fftw_destroy_plan(pf_y);
        fftw_destroy_plan(pf_s);
        fftw_destroy_plan(pb_x);
        fftw_destroy_plan(pb_y);
        fftw_destroy_plan(pb_s);
        delete []ux;
        delete []uy;
        delete []s;
        delete []uxf;
        delete []uyf;
        delete []sf;
    }
    
    std::cout<<"U After Diff, top:"<<U1[Idx(N[0]/3-10,N[1]/2+10)]<<",btm:"<<U1[2*Idx(N[0]/3+10,N[1]/2)]<<std::endl;
}

template<class T, int n>
Vec<T,n> abs(Vec<T,n> v) {
    for(int i=0; i<n; i++) v[i] = std::abs(v[i]);
    return v;
}

template<class T, int n>
void Fluid<T,n>::Project(){
    // U1, U0, dt
    if(!FFT_scheme){
        T dx2 = D[0]*D[0], dy2 = D[1]*D[1], d2 = dx2*dy2;
        
        // divU
        for(int i=0; i<N[0]; i++){
            for(int j=0; j<N[1]; j++){
                div[Idx(i,j)] = 0.5*((U1[Idx(i+1,j)][0]-U1[Idx(i-1,j)][0])/D[0] + (U1[Idx(i,j+1)][1]-U1[Idx(i,j-1)][1])/D[1]);
                P[Idx(i,j)] = 0.0;
            }
        }
        
        // vorticity
        std::vector<T> uCurls(N[0]*N[1]); 
        for(int i=0; i<N[0]; i++){
            for(int j=0; j<N[1]; j++){
                uCurls[Idx(i,j)] = 0.5*((U1[Idx(i+1,j)][1]-U1[Idx(i-1,j)][1])/D[0] + (-U1[Idx(i,j+1)][0]+U1[Idx(i,j-1)][0])/D[1]);
            }
        }
        for(int i=0; i<N[0]; i++) {
            vec mx = {1e3, 1e3}, mn = {-1e3, -1e3};
            for(int j=0; j<N[1]; j++) {
                T cl = uCurls[Idx(i-1, j)];
                T cr = uCurls[Idx(i+1, j)];
                T cb = uCurls[Idx(i, j-1)];
                T ct = uCurls[Idx(i, j+1)];
                T cc = uCurls[Idx(i, j)];
                vec force = vec(abs(ct) - abs(cb), abs(cl) - abs(cr));
                force = 1e-2 * cc * 7.0 * force.normalized(); // 7.0 factor
                U1[Idx(i,j)] = componentwise_min(componentwise_max(U1[Idx(i,j)] + force * dt, mn), mx);
            }
        }

        // calculate P
        int i0 = N[0]/2, j0 = N[1]/2;
        T prev = 0.0;
        for(int iter=0; iter<20; iter++){
            prev = P[Idx(i0,j0)];
            for(int i=0; i<N[0]; i++){
                for(int j=0; j<N[0]; j++){
                    T A = (P[Idx(i+1,j)]+P[Idx(i-1,j)])*dy2;
                    T B = (P[Idx(i,j+1)]+P[Idx(i,j-1)])*dx2;
                    P[Idx(i,j)] = (A+B-div[Idx(i,j)]*d2)/(2.0*(dx2+dy2));
                }
            }
            if(iter==19 && P[Idx(i0,j0)]-prev>small_t) std::cout<<"Not converge: "<<P[Idx(i0,j0)]-prev<<std::endl;
        }

        // w4 = w3 - divP
        for(int i=0; i<N[0]; i++){
            for(int j=0; j<N[1]; j++){
                vec gradP(0.5*(P[Idx(i+1,j)]-P[Idx(i-1,j)])/D[0], 0.5*(P[Idx(i,j+1)]-P[Idx(i,j-1)])/D[1]);
                U1[Idx(i,j)] -= gradP;
            }
        }
        std::cout<<"U After P, top:"<<U1[Idx(N[0]/4,N[1]/2)]<<",btm:"<<U1[Idx(3*N[0]/4,N[1]/2)]<<std::endl;
    }
    // source: disspation term
    T cs = 1.0 + dt*aS;
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            S1[Idx(i,j)] = S1[Idx(i,j)]/cs;
        }
    }
    std::cout<<"S After P, top:"<<S1[Idx(N[0]/4,N[1]/2)]<<",btm:"<<S1[Idx(3*N[0]/4,N[1]/2)]<<std::endl;
}

// RK2
template<class T, int n>
inline void Fluid<T,n>::TraceParticle(vec& X1, vec& X0)
{
    vec k1(Interpolate(U0,X1));
    vec Xt = X1 - dt*k1;
    
    vec k2(Interpolate(U0,Xt));
    X0 = X1 - (k1+k2)/2*dt;
}


template<class T, int n>
inline Vec<T,n> Fluid<T,n>::Interpolate(std::vector<vec>& U, vec& X){
    // portion along axis
    T x = N[1]/2+X[0]/L[1]*N[1], y = N[0]/2-X[1]/L[0]*N[0];
    
    int i0 = y<N[0]/2? floor(y): floor(y)-1, i1 = i0 + 1,
        j0 = x<N[1]/2? floor(x): floor(x)-1, j1 = j0 + 1;

    // x = (x-0.5-j0), y = (y-0.5-i0);
    x = std::min(std::max((x-0.5-j0),0.0), 1.0);
    y = std::min(std::max((y-0.5-i0),0.0), 1.0);
    
    vec ut = (1-x)*U[Idx(i0,j0)] + x*U[Idx(i0,j1)],
        ub = (1-x)*U[Idx(i1,j0)] + x*U[Idx(i1,j1)];
    
    return (1-y)*ut + y*ub;
}


template<class T, int n>
inline T Fluid<T,n>::Interpolate(std::vector<T>& S, vec& X){
    // portion along axis
    T x = N[1]/2+X[0]/L[1]*N[1], y = N[0]/2-X[1]/L[0]*N[0];
    
    int i0 = y<N[0]/2? floor(y): floor(y)-1, i1 = i0 + 1,
        j0 = x<N[1]/2? floor(x): floor(x)-1, j1 = j0 + 1;

    x = std::min(std::max((x-0.5-j0),0.0), 1.0);
    y = std::min(std::max((y-0.5-i0),0.0), 1.0);

    T st = (1-x)*S[Idx(i0,j0)] + x*S[Idx(i0,j1)],
      sb = (1-x)*S[Idx(i1,j0)] + x*S[Idx(i1,j1)];

    return (1-y)*st + y*sb;
}

template<class T, int n>
inline int Fluid<T,n>::XtoIdx(const vec& X){
    int j = N[1]/2+X[0]/L[1]*N[1], i = N[0]/2-X[1]/L[0]*N[0];
    return Idx(i,j);
}

template<class T, int n>
inline T Fluid<T,n>::IdxtoX(const int i, const int dim){
    if(dim==0) return (0.5 - (i+0.5)/N[dim])*L[dim];
    else return ((i+0.5)/N[dim] - 0.5)*L[dim];
}

template<class T, int n>
inline int Fluid<T,n>::Idx(int i, int j){
    i = (i+N[0]) % N[0];
    j = (j+N[1]) % N[1];
    return i*N[1] + j;
}

template<class T, int n>
void Fluid<T,n>::display()
{
    std::cout<<"display"<<std::endl;
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            int val = std::min(S0[Idx(i,j)], (T)255);
            image_color[Idx(i,j)] = make_pixel(val, val, val);
        }
    }
    // obstacle
    if(obstacle){
        for(int i=0; i<N[0]; i++)
            for(int j=0; j<N[1]; j++)
                if((i-N[0]/2)*(i-N[0]/2)+(j-N[1]/2)*(j-N[1]/2)<std::min(N[0],N[1]))
                    image_color[Idx(i,j)] = make_pixel(200, 133, 20);
    }
    dump_png(image_color,N[1],N[0],"output.png");
}

template<class T, int n>
void Fluid<T,n>::AddSource(const char* filename){
    Read_png(image_color,N[1],N[0],filename);

    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            int r, g, b;
            from_pixel(image_color[Idx(i,j)],r,g,b);
            S0[Idx(i,j)] += (r+g+b)/3 * 0.9;
        }
    }
}


template<class T, int n>
void Fluid<T,n>::Check_Symmetry(std::vector<vec>& U, std::vector<T>& S, 
    const std::string& s){
    if(!U.empty()){
        for(int i=0; i<N[0]; i++){
            for(int j=0; j<N[1]; j++){
                if(U[Idx(i,j)][0]!=U[Idx(N[0]-1-i,j)][0 || 
                    U[Idx(i,j)][1]!=U[Idx(N[0]-1-i,j)][1]]){
                    std::cout<<"U h:"<<s<<std::endl; break;
                }
                if(U[Idx(i,j)][0]!=U[Idx(i,N[1]-1-j)][0] || 
                    U[Idx(i,j)][1]!=U[Idx(i,N[1]-1-j)][1]){
                    std::cout<<"U v:"<<s<<std::endl; break;
                }
            }
        }
    }
    if(!S.empty()){
        for(int i=0; i<N[0]; i++){
            for(int j=0; j<N[1]; j++){
                if(S[Idx(i,j)]!=S[Idx(N[0]-1-i,j)]){
                    std::cout<<"S h:"<<s<<" row:"<<i<<std::endl; break;
                }
                if(S[Idx(i,j)]!=S[Idx(i,N[1]-1-j)]){
                    std::cout<<"S v:"<<s<<std::endl; break;
                }
            }
        }
    }
}
