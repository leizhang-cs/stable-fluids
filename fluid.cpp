#include <fftw3.h>
#include <cassert>
#include "fluid.h"
#include "pixel.h"


// I would test uniform translation (easy to debug) and the Taylor-Green vortex artefacts near the edges.
using Type = double;
const int Ndim = 2;

template class Fluid<Type,Ndim>;
using vec = Vec<Type,Ndim>;
using pixel = unsigned int;

static const double small_t = 1e-4;


template<class T, int n>
void Fluid<T,n>::simulate(vec& F, T Source, vec& X){
    Vstep(F, Source, X);
    //S.step();
    std::swap(U0, U1);
    std::swap(S0, S1);
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

    //int index = XtoIdx(X);

    // w1 = f(w0)
    vec du_top = F*dt;
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            U1[Idx(i,j)] = du_top;
            //U1[Idx(i,j)] = F*dt*S0[Idx(i,j)];
        }
    }
    
    vec du_btm = F*dt;
    for(int i=N[0]/2; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            U1[Idx(i,j)] = du_btm;
            //U1[Idx(i,j)] = F*dt*S0[Idx(i,j)];
        }
    }
    
    // add source
    T ds = S;
    for(int i=3*N[0]/8; i<5*N[0]/8; i++){
        for(int j=0*N[1]/8; j<2*N[1]/8; j++){
            S1[Idx(i,j)] = ds;
        }
    }
    // obstacle
    if(obstacle){
        for(int i=0; i<N[0]; i++)
            for(int j=0; j<N[1]; j++)
                if((i-64)*(i-64)+(j-64)*(j-64)<81){
                    U0[Idx(i,j)] = vec(0,0);         
                }
    }
    //Check_Symmetry(U1, S1, "addF");
    std::cout<<"du top: "<<du_top<<", btm: "<<du_btm<<std::endl;
    std::cout<<"delta source: "<<ds<<std::endl;
}

template<class T, int n>
inline T Fluid<T,n>::IdxtoX(const int i, const int dim){
    if(dim==0) return (0.5 - (i+0.5)/N[dim])*L[dim];
    else return ((i+0.5)/N[dim] - 0.5)*L[dim];
}

template<class T, int n>
void Fluid<T,n>::Advect(){
    // U1, U0, dt
    // TraceParticle: method of charactristic
    std::vector<vec> Ut(N[0]*N[1]);

    for(int i=0; i<N[0]; i++){
        T y = IdxtoX(i,0); // (0, N[1]) <=> (L[1], -L[1])
        //std::cout<<i<<"=>"<<y<<"\n";
        for(int j=0; j<N[1]; j++){
            vec X0;
            vec X1(IdxtoX(j,1), y); // (0, N[0]) <=> (-L[0], L[0])
            TraceParticle(X1, X0);
            //if(i%17==0) std::cout<<j<<"=>"<<X1[0]<<"\n";
            //if(j==17) std::cout<<X1<<"=>"<<X0<<"\n";
            // w2 = f(w1)
            // Interpolation
            vec ut;
            T st;
            int flag = 2;
            //if((i==5*N[1]/8+10) && j==4*N[1]/8){ std::cout<<X1<<"=>"<<X0<<"\n"; std::cout<<"i,j: "<<i<<" "<<j<<" "; debug_flag = true; }
            Interpolate(U0, ut, S0, st, X0, flag);
            debug_flag = false;
            
            U1[Idx(i,j)] += ut;
            S1[Idx(i,j)] += st;
        }
    }
    //Check_Symmetry(U1, S1, "adv");
    std::cout<<"U After Adv, top:"<<U1[Idx(N[0]/2-10,N[1]/2)]<<",btm:"<<U1[Idx(N[0]/2+10,N[1]/2)]<<std::endl;
    //std::cout<<"S after Adv: "<<min_val<<std::endl;
    //std::cout<<"Sssssssss:"<<S1[Idx(100,64)]<<"\n";
}


template<class T, int n>
void Fluid<T,n>::Diffuse(){
    // U1, U0, visc, dt 
    // conjugate gradient. FTCS. BTCS???
    // unknown: U1[Idx(i,j)]
    // FFT
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
    std::cout<<"U before Diff, top:"<<U1[Idx(N[0]/2-10,N[1]/2)]<<",btm:"<<U1[Idx(N[0]/2+10,N[1]/2)]<<std::endl;
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
            //T f = 1/(1 + visc*dt*r);
            T f = exp(-visc*dt*r);
            if(j%29==0) std::cout<<f<<std::endl;
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
    std::cout<<"U After Proj, top:"<<U1[Idx(N[0]/2-10,N[1]/2)]<<",btm:"<<U1[Idx(N[0]/2+10,N[1]/2)]<<std::endl;
    //std::cout<<"S after Diff: "<<min_val<<std::endl;
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
    /*
    // numerical method
    T k = visc*dt; // sign -1?
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            // TODO: optimization
            vec ux = (U0[Idx(i+1,j)] - 2.0*U0[Idx(i,j)] + U0[Idx(i-1,j)])/(D[0]*D[0]);
            vec uy = (U0[Idx(i,j+1)] - 2.0*U0[Idx(i,j)] + U0[Idx(i,j-1)])/(D[1]*D[1]);
            U1[Idx(i,j)] += k*(ux + uy);
        }
    }
    if(!pbc) boundry_condition(U1);
    */
}

template<class T, int n>
void Fluid<T,n>::Project(){
    // U1, U0, dt
    
    T dx2 = D[0]*D[0], dy2 = D[1]*D[1], d2 = dx2*dy2;
    
    // divU
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            div[Idx(i,j)] = 0.5*((U1[Idx(i+1,j)][0]-U1[Idx(i-1,j)][0])/D[0] + (U1[Idx(i,j+1)][1]-U1[Idx(i,j-1)][1])/D[1]);
            P[Idx(i,j)] = 0.0;
        }
    }
    if(!pbc) boundry_condition(div);

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
        if(!pbc) boundry_condition(P);
        if(iter==19 && P[Idx(i0,j0)]-prev>small_t) std::cout<<"Not converge: "<<P[Idx(i0,j0)]-prev<<std::endl;
    }

    // w4 = w3 - divP
    T max_U = 0;
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            vec gradP(0.5*(P[Idx(i+1,j)]-P[Idx(i-1,j)])/D[0], 0.5*(P[Idx(i,j+1)]-P[Idx(i,j-1)])/D[1]);
            U1[Idx(i,j)] -= gradP;
            max_U = std::max(max_U, U1[Idx(i,j)][0]);
        }
    }
    if(!pbc) boundry_condition(U1);
    std::cout<<"U After P, top:"<<U1[Idx(N[0]/4,N[1]/2)]<<",btm:"<<U1[Idx(3*N[0]/4,N[1]/2)]<<std::endl;
    

    // source: disspation term
    T max_val = 0;
    T cs = 1.0 + dt*aS;
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            S1[Idx(i,j)] = S1[Idx(i,j)]/cs;
        }
    }
    
    // std::cout<<"S after P:"<<max_val<<std::endl;
}

template<class T, int n>
inline void Fluid<T,n>::TraceParticle(const vec& X1, vec& X0)
{
    vec k1(Interpolate(U0,X0));
    vec Xt = X1 - dt*k1;
    vec k2(Interpolate(U0,Xt));
    X0 = X1 - (k1+k2)/2*dt;

    X0 = Xt;
}

// TODO: write as mode 1 and 2
template<class T, int n>
inline void Fluid<T,n>::Interpolate(const std::vector<vec>& U, vec& ret_U, 
    const std::vector<T>& S, T& ret_S, vec& X, int flag){
    // portion along axis
    // T x = X[0]/L[1], y = X[1]/L[0];
    // int i0 = floor(N[0]/2-y*N[0]), i1 = i0<N[0]/2? i0+1: i0-1,
    //     j0 = floor(N[1]/2+x*N[1]), j1 = j0<N[1]/2? j0+1: j0-1;
    // if(i0>i1) std::swap(i0,i1);
    // if(j0>j1) std::swap(j0,j1);
    // // portion in a cell
    // x = (N[1]/2+x*N[1]-j0); y = (N[0]/2-y*N[0]-i0);
    // //x = 0.5, y = 0.5;
    // if(x<0 || x>1) std::cout<<j0<<" "<<N[1]/2+ X[0]/L[1]*N[1]<<std::endl;
    T x = N[1]/2+X[0]/L[1]*N[1], y = N[0]/2-X[1]/L[0]*N[0];
    
    int i0 = y<N[0]/2? floor(y): floor(y)-1, i1 = i0 + 1,
        j0 = x<N[1]/2? floor(x): floor(x)-1, j1 = j0 + 1;

    if(debug_flag) std::cout<<y<<" "<<i0<<" "<<i1<<"  "<<x<<"  "<<j0<<" "<<j1<<"\n";
    if(debug_flag) std::cout<<x-0.5<<"  "<<y-0.5<<std::endl;

    x = (x-0.5-j0), y = (y-0.5-i0);
    
    if(debug_flag) std::cout<<"portion: "<<x<<"  "<<y<<std::endl;
    if(debug_flag) std::cout<<S[Idx(i0,j0)]<<" "<<S[Idx(i1,j0)]<<" "<<S[Idx(i0,j1)]<<" "<<S[Idx(i1,j1)]<<"\n";

    if(flag==0 || flag==2){
        vec ut = (1-x)*U[Idx(i0,j0)] + x*U[Idx(i0,j1)],
            ub = (1-x)*U[Idx(i1,j0)] + x*U[Idx(i1,j1)];
        ret_U = (1-y)*ut + y*ub;
    }
    if(flag==1 || flag==2){
        T st = (1-x)*S[Idx(i0,j0)] + x*S[Idx(i0,j1)],
          sb = (1-x)*S[Idx(i1,j0)] + x*S[Idx(i1,j1)];
        ret_S = (1-y)*st + y*sb;
        if(debug_flag) std::cout<<j0<<" "<<j1<<"  S:"<<st<<"\n";
    }
    // if(debug_flag){
    //     x = X[0]/L[0];
    //     std::cout<<X<<" S"<<ret_S<<" i,j:"<<j0<<" "<<j1<<"   "<<N[1]/2-x*N[1]<<" floor:"<<floor(N[0]/2-x*N[0])<<"\n";
    // }
}

template<class T, int n>
inline Vec<T,n> Fluid<T,n>::Interpolate(std::vector<vec>& U, vec& X){
    // portion along axis
    // T x = X[0]/L[1]+0.5, y = X[1]/L[0]+0.5;
    // int i0 = floor(y*N[0]), i1 = (i0+1)%N[0],
    //     j0 = floor(x*N[1]), j1 = (j0+1)%N[1];
    // // portion in a cell
    // x = (x*N[1]-j0); y = (y*N[0]-i0);
    
    // vec ut = (1-x)*U[Idx(i0,j0)] + x*U[Idx(i0,j1)],
    //     ub = (1-x)*U[Idx(i1,j0)] + x*U[Idx(i1,j1)];
    // return vec((1-y)*ut + y*ub);
    T x = N[1]/2+X[0]/L[1]*N[1], y = N[0]/2-X[1]/L[0]*N[0];
    
    int i0 = y<N[0]/2? floor(y): floor(y)-1, i1 = i0 + 1,
        j0 = x<N[1]/2? floor(x): floor(x)-1, j1 = j0 + 1;

    x = (x-0.5-j0), y = (y-0.5-i0);  
    
    vec ut = (1-x)*U[Idx(i0,j0)] + x*U[Idx(i0,j1)],
        ub = (1-x)*U[Idx(i1,j0)] + x*U[Idx(i1,j1)];
    return (1-y)*ut + y*ub;
    
}


template<class T, int n>
inline T Fluid<T,n>::Interpolate(std::vector<T>& S, vec& X){
    // portion along axis
    T x = X[0]/L[1]+0.5, y = X[1]/L[0]+0.5;
    int i0 = floor(y*N[0]), i1 = (i0+1)%N[0],
        j0 = floor(x*N[1]), j1 = (j0+1)%N[1];
    // portion in a cell
    x = (x*N[1]-j0); y = (y*N[0]-i0);
    T st = (1-x)*S[Idx(i0,j0)] + x*S[Idx(i0,j1)],
      sb = (1-x)*S[Idx(i1,j0)] + x*S[Idx(i1,j1)];
    return (1-y)*st + y*sb;
}

template<class T, int n>
inline int Fluid<T,n>::XtoIdx(const vec& X){
    int j = (0.5-X[0]/L[1])*N[1] ,
        i = (0.5-X[1]/L[0])*N[0] ;
    return Idx(i,j);
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
            //if(S0[Idx(i,j)]<0){ std::cout<<"<0"<<"\n"; }
            int val = std::min(S0[Idx(i,j)], (T)255);
            image_color[Idx(i,j)] = make_pixel(val, val, val);
        }
    }
    // obstacle
    if(obstacle){
        for(int i=0; i<N[0]; i++)
            for(int j=0; j<N[1]; j++)
                if((i-64)*(i-64)+(j-64)*(j-64)<81)
                    image_color[Idx(i,j)] = make_pixel(200, 133, 20);
    }
    dump_png(image_color,N[0],N[1],"output.png");
}

template<class T, int n>
void Fluid<T,n>::AddSource(const char* filename){
    //TODO readpng modify out
    Read_png(image_color,N[0],N[1],filename);

    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            int r, g, b;
            from_pixel(image_color[Idx(i,j)],r,g,b);
            S0[Idx(i,j)] += (r+g+b)/3 * 0.5;
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

template<class T, int n>
void Fluid<T,n>::boundry_condition(std::vector<T>& var){
    var[Idx(0,0)] = var[Idx(1,1)];
    var[Idx(0,1)] = var[Idx(1,1)];
    var[Idx(1,0)] = var[Idx(1,1)];
    var[Idx(N[0]-1,N[0]-1)] = var[Idx(N[0]-2,N[0]-2)];
    var[Idx(N[0]-1,N[0]-2)] = var[Idx(N[0]-2,N[0]-2)];
    var[Idx(N[0]-2,N[0]-1)] = var[Idx(N[0]-2,N[0]-2)];
}

template<class T, int n>
void Fluid<T,n>::boundry_condition(std::vector<vec>& var){
    var[Idx(0,0)] = var[Idx(1,1)];
    var[Idx(0,1)] = var[Idx(1,1)];
    var[Idx(1,0)] = var[Idx(1,1)];
    var[Idx(N[0]-1,N[0]-1)] = var[Idx(N[0]-2,N[0]-2)];
    var[Idx(N[0]-1,N[0]-2)] = var[Idx(N[0]-2,N[0]-2)];
    var[Idx(N[0]-2,N[0]-1)] = var[Idx(N[0]-2,N[0]-2)];
}

