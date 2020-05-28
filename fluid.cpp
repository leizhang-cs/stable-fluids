#include "fluid.h"
#include "pixel.h"
#include <fftw3.h>


template class Fluid<double,2>;
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
void Fluid<T,n>::display()
{
    std::cout<<"display"<<std::endl;
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            int val = std::min(S0[Idx(i,j)], (T)255);
            image_color[Idx(i,j)] = make_pixel(val, val, val);
        }
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
            S0[Idx(i,j)] = (r+g+b)/3 * 0.5;
        }
    }
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
    vec du_top = F*dt/(density*L[0]*L[1]);
    for(int i=0; i<N[0]/2; i++){
        for(int j=0; j<N[1]; j++){
            U1[Idx(i,j)] = du_top;
        }
    }
    
    F[0] = -F[0];
    vec du_btm = F*dt/(density*L[0]*L[1]);
    for(int i=N[0]/2; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            U1[Idx(i,j)] = du_btm;
        }
    }
    
    // add source
    T ds = S/(N[0]*N[1])*16;
    for(int i=3*N[0]/8; i<5*N[0]/8; i++){
        for(int j=0; j<N[1]/4; j++){
            S1[Idx(i,j)] = ds;
        }
    }
    
    std::cout<<"du top: "<<du_top<<", btm: "<<du_btm<<std::endl;
    std::cout<<"delta source: "<<ds<<std::endl;
}

template<class T, int n>
void Fluid<T,n>::Advect(){
    // U1, U0, dt
    // TraceParticle: method of charactristic
    std::vector<vec> Ut(N[0]*N[1]);

    T max_val = 0;
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            vec X1, X0;
            vec index(i, j);
            X1 = O + index*D;
            TraceParticle(X1, X0);
            
            // w2 = f(w1)
            U1[Idx(i,j)] += U0[XtoIdx(X0)];
            S1[Idx(i,j)] += S0[XtoIdx(X0)];
            max_val = std::max(max_val, S1[Idx(i,j)]);
        }
    }

    /*for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            U1[Idx(i,j)] = Ut[Idx(i,j)];
            S1[Idx(i,j)] = St[Idx(i,j)];
        }
    }*/
    //std::cout<<"U After Adv, top:"<<U1[Idx(N[0]/4,N[1]/2)]<<",btm:"<<U1[Idx(3*N[0]/4,N[1]/2)]<<std::endl;
    //std::cout<<"S After Adv: "<<max_val<<std::endl;
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

    T max_val = 0;
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            ux[Idx(i,j)][0] = U1[Idx(i,j)][0];
            uy[Idx(i,j)][0] = U1[Idx(i,j)][1];
            s[Idx(i,j)][0] = S1[Idx(i,j)];
            max_val = std::max(max_val, S1[Idx(i,j)]);
        }
    }
    std::cout<<"U before Diff, top:"<<U1[Idx(N[0]/4,N[1]/2)]<<",btm:"<<U1[Idx(3*N[0]/4,N[1]/2)]<<std::endl;
    //std::cout<<"S before Diff: "<<max_val<<std::endl;
    fftw_plan pf_x = fftw_plan_dft_1d(N[0]*N[1], ux, uxf, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan pf_y = fftw_plan_dft_1d(N[0]*N[1], uy, uyf, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan pf_s = fftw_plan_dft_1d(N[0]*N[1], s, sf, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(pf_x);
    fftw_execute(pf_y);
    fftw_execute(pf_s);

    // wave length
    T lambda1 = L[0], lambda2 = L[1];
    // wave number
    // k(k1,k2), k1 = 2*pi/lambda1, k2 = 2*pi/lambda2
    vec k(2*pi/lambda1,2*pi/lambda2);
    T k_2 = k.magnitude_squared();
    T c = 1.0 + visc*dt*k_2;
    T cs = 1.0 + kS*dt*k_2;
    std::vector<vec> ux_P(num), uy_P(num);
    //std::cout<<"coeff:"<<c<<std::endl;

    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            ux_P[Idx(i,j)][0] = uxf[Idx(i,j)][0] /= c;
            ux_P[Idx(i,j)][1] = uxf[Idx(i,j)][1] /= c;
            uy_P[Idx(i,j)][0] = uyf[Idx(i,j)][0] /= c;
            uy_P[Idx(i,j)][1] = uyf[Idx(i,j)][1] /= c;
            sf[0][Idx(i,j)] /= cs;
            sf[1][Idx(i,j)] /= cs;
        }
    }
    /*
    // Projection by FFT
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            vec v;
            v = (k[0]*k[0]*ux_P[Idx(i,j)] + k[0]*k[1]*(ux_P[Idx(i,j)]+uy_P[Idx(i,j)]) +
                k[1]*k[1]*uy_P[Idx(i,j)]) / k_2;
            uxf[Idx(i,j)][0] -= v[0]; uxf[Idx(i,j)][1] -= v[1];
            uyf[Idx(i,j)][0] -= v[0]; uyf[Idx(i,j)][1] -= v[1];
        }
    }
    */

    fftw_plan pb_x = fftw_plan_dft_1d(N[0]*N[1], uxf, ux, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan pb_y = fftw_plan_dft_1d(N[0]*N[1], uyf, uy, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan pb_s = fftw_plan_dft_1d(N[0]*N[1], sf, s, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(pb_x);
    fftw_execute(pb_y);
    fftw_execute(pb_s);

    max_val = 0;
    T max_U = 0;
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            U1[Idx(i,j)][0] = ux[Idx(i,j)][0]/num;
            U1[Idx(i,j)][1] = uy[Idx(i,j)][0]/num;
            S1[Idx(i,j)] = s[Idx(i,j)][0]/num;
            max_val = std::max(max_val, S1[Idx(i,j)]);
            max_U = std::max(max_U, U1[Idx(i,j)][0]);
        }
    }    
    std::cout<<"U After Diff, top:"<<U1[Idx(N[0]/4,N[1]/2)]<<",btm:"<<U1[Idx(3*N[0]/4,N[1]/2)]<<std::endl;
    //std::cout<<"S after Diff: "<<max_val<<std::endl;

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
    /*
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
    std::cout<<"U after P:"<<max_U<<std::endl;
    */

    // source: disspation term
    T max_val = 0;
    T cs = 1.0 + dt*aS;
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            S1[Idx(i,j)] = S1[Idx(i,j)]/cs;
            max_val = std::max(max_val, S1[Idx(i,j)]);
        }
    }
    
    std::cout<<"S after P:"<<max_val<<std::endl;
}

template<class T, int n>
inline void Fluid<T,n>::TraceParticle(vec& X1, vec& X0)
{
    vec k1 = U0[XtoIdx(X1)];
    vec X_t = X1 - dt*k1;
    // TODO: t-dt??
    vec k2 = U0[XtoIdx(X_t)];
    X0 = X1 - (k1+k2)/2*dt;
}


template<class T, int n>
inline int Fluid<T,n>::XtoIdx(vec& X){
    int i = (X[0]+L[0])*N[0]/L[0] + 0.5,
        j = (X[1]+L[1])*N[1]/L[1] + 0.5;
    return Idx(i,j);
}

template<class T, int n>
inline int Fluid<T,n>::Idx(int i, int j){
    i = (i+N[0]) % N[0];
    j = (j+N[1]) % N[1];
    return i*N[1] + j;
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



/*
inline pixel make_pixel(int r, int g, int b)
{
    return (r<<24)|(g<<16)|(b<<8)|0xff;
}

// Dump an image to file.
void dump_png(pixel* data,int width,int height,const char* filename)
{
    FILE* file=fopen(filename,"wb");
    assert(file);

    png_structp png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
    assert(png_ptr);
    png_infop info_ptr=png_create_info_struct(png_ptr);
    assert(info_ptr);
    bool result=setjmp(png_jmpbuf(png_ptr));
    assert(!result);
    png_init_io(png_ptr,file);
    int color_type=PNG_COLOR_TYPE_RGBA;
    png_set_IHDR(png_ptr,info_ptr,width,height,8,color_type,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);

    pixel** row_pointers=new pixel*[height];
    for(int j=0;j<height;j++) row_pointers[j]=data+width*(height-j-1);
    png_set_rows(png_ptr,info_ptr,(png_byte**)row_pointers);
    png_write_png(png_ptr,info_ptr,PNG_TRANSFORM_BGR|PNG_TRANSFORM_SWAP_ALPHA,0);
    delete[] row_pointers;
    png_destroy_write_struct(&png_ptr,&info_ptr);
    fclose(file);
}*/