#include "fluid.h"
#include "pixel.h"
#include <fftw3.h>


template class Fluid<double,2>;
using pixel = unsigned int;

static double small_t = 1e-4;

template<class T, int n>
void Fluid<T,n>::simulate(vec& F, T Source, vec& X){
    Vstep(F,X);
    //S.step();
    std::swap(U0, U1);
    //std::swap(S0, S1);
}

template<class T, int n>
void Fluid<T,n>::display()
{
    std::cout<<"display"<<std::endl;
    for(int i=0; i<size[0]; i++){
        for(int j=0; j<size[1]; j++){
            int val = std::min(U0[Idx(i,j)].magnitude(), (T)255);
            image_color[Idx(i,j)] = make_pixel(val, val, val);
        }
    }
    
    dump_png(image_color,size[0],size[1],"output.png");
}

template<class T, int n>
void Fluid<T,n>::Vstep(vec& F, vec& X){
    // U0, U1, visc, F, dt
    AddForce(F, X);
    Advect();
    Diffuse();
    Project();
}

template<class T, int n>
void Fluid<T,n>::AddForce(vec F, vec X){
    // U1 = U0 + F*dt
    // TODO: volume. dm: delta momentum
    vec dm = F*dt/density;

    int index = XtoIdx(X);
    U1[index] = U0[index] + dm;
    std::cout<<U1[index]<<std::endl;
}

template<class T, int n>
void Fluid<T,n>::Advect(){
    // U1, U0, dt
    // TraceParticle: method of charactristic
    for(int i=0; i<size[0]; i++){
        for(int j=0; j<size[1]; j++){
            vec X1, X0;
            vec index1(i, j);
            X1 = O + index1*D;
            TraceParticle(X1, X0);
            // interpolate
            U1[Idx(i,j)] += U0[XtoIdx(X0)];
        }
    }
    
}

template<class T, int n>
void Fluid<T,n>::Diffuse(){
    // U1, U0, visc, dt 
    // conjugate gradient. FTCS. BTCS???
    // unknown: U1[Idx(i,j)]
    T k = visc*dt; // sign -1?
    // U1[Idx(0,0)] += k*(U0[Idx(1,0)] - U0[Idx(0,0)])/(D[0]*D[0]) + (U0[Idx(0,1)] - U0[Idx(0,0)])/(D[1]*D[1]);
    // U1[Idx(1,0)] += k*(U0[2][0] - 2.0*U0[Idx(1,0)] + U0[Idx(0,0)])/(D[0]*D[0]) + (U0[Idx(1,1)] - U0[Idx(1,0)])/(D[1]*D[1]);
    // U1[Idx(0,1)] += k*(U0[Idx(1,1)] - U0[Idx(0,1)])/(D[0]*D[0]) + (U0[0][2] - 2.0*U0[Idx(0,1)] + U0[Idx(0,0)])/(D[1]*D[1]);
    for(int i=0; i<size[0]; i++){
        for(int j=0; j<size[1]; j++){
            // TODO: optimization
            vec ux = (U0[Idx(i+1,j)] - 2.0*U0[Idx(i,j)] + U0[Idx(i-1,j)])/(D[0]*D[0]);
            vec uy = (U0[Idx(i,j+1)] - 2.0*U0[Idx(i,j)] + U0[Idx(i,j-1)])/(D[1]*D[1]);
            U1[Idx(i,j)] += k*(ux + uy);
        }
    }
    if(!pbc) boundry_condition(U1);
}

template<class T, int n>
void Fluid<T,n>::Project(){
    // U1, U0, dt
    T dx2 = D[0]*D[0], dy2 = D[1]*D[1], d2 = dx2*dy2;
    
    // divU
    for(int i=0; i<size[0]; i++){
        for(int j=0; j<size[1]; j++){
            div[Idx(i,j)] = 0.5*((U1[Idx(i+1,j)][0]-U1[Idx(i-1,j)][0])/D[0] + (U1[Idx(i,j+1)][1]-U1[Idx(i,j-1)][1])/D[1]);
            P[Idx(i,j)] = 0.0;
        }
    }
    if(!pbc) boundry_condition(div);

    // calculate P
    int i0 = size[0]/2, j0 = size[1]/2;
    T prev = 0.0;
    for(int iter=0; iter<20; iter++){
        prev = P[Idx(i0,j0)];
        for(int i=0; i<size[0]; i++){
            for(int j=0; j<size[0]; j++){
                T A = (P[Idx(i+1,j)]+P[Idx(i-1,j)])*dy2;
                T B = (P[Idx(i,j+1)]+P[Idx(i,j-1)])*dx2;
                P[Idx(i,j)] = (A+B-div[Idx(i,j)]*d2)/(2.0*(dx2+dy2));
            }
        }
        if(!pbc) boundry_condition(P);
        if(iter==19 && P[Idx(i0,j0)]-prev>small_t) std::cout<<"Not converge: "<<P[Idx(i0,j0)]-prev<<std::endl;
    }

    // w4 = w3 - divP
    T max_val = 0;
    for(int i=0; i<size[0]; i++){
        for(int j=0; j<size[1]; j++){
            vec gradP(0.5*(P[Idx(i+1,j)]-P[Idx(i-1,j)])/D[0], 0.5*(P[Idx(i,j+1)]-P[Idx(i,j-1)])/D[1]);
            U1[Idx(i,j)] -= gradP;
            max_val = std::max(max_val, U1[Idx(i,j)].magnitude());
        }
    }
    std::cout<<"U after P:"<<max_val<<std::endl;
    if(!pbc) boundry_condition(U1);
}

template<class T, int n>
void Fluid<T,n>::TraceParticle(vec& X1, vec& X0)
{
    // TODO RK
    X0 = X1 - dt*U0[XtoIdx(X1)];
}

template<class T, int n>
Vec<T,n> Fluid<T,n>::Interpolate(std::vector<vec>& U, vec& X)
{
    // X U
    int i = std::min(static_cast<int>(X[0]*L[0]/D[0]), size[0]-1), 
        j = std::min(static_cast<int>(X[1]*L[1]/D[1]), size[1]-1);
    return U[Idx(i,j)];
}

template<class T, int n>
int Fluid<T,n>::Idx(int i, int j){
    i = (i+size[0])%size[0];
    j = (j+size[1])%size[1];
    return i*size[1] + j;
}


template<class T, int n>
int Fluid<T,n>::XtoIdx(vec& X){
    int i = std::min(static_cast<int>(X[0]*L[0]/D[0]), size[0]-1),
        j = std::min(static_cast<int>(X[1]*L[1]/D[1]), size[1]-1);
    return Idx(i,j);
}

template<class T, int n>
void Fluid<T,n>::boundry_condition(std::vector<T>& var){
    var[Idx(0,0)] = var[Idx(1,1)];
    var[Idx(0,1)] = var[Idx(1,1)];
    var[Idx(1,0)] = var[Idx(1,1)];
    var[Idx(size[0]-1,size[0]-1)] = var[Idx(size[0]-2,size[0]-2)];
    var[Idx(size[0]-1,size[0]-2)] = var[Idx(size[0]-2,size[0]-2)];
    var[Idx(size[0]-2,size[0]-1)] = var[Idx(size[0]-2,size[0]-2)];
}

template<class T, int n>
void Fluid<T,n>::boundry_condition(std::vector<vec>& var){
    var[Idx(0,0)] = var[Idx(1,1)];
    var[Idx(0,1)] = var[Idx(1,1)];
    var[Idx(1,0)] = var[Idx(1,1)];
    var[Idx(size[0]-1,size[0]-1)] = var[Idx(size[0]-2,size[0]-2)];
    var[Idx(size[0]-1,size[0]-2)] = var[Idx(size[0]-2,size[0]-2)];
    var[Idx(size[0]-2,size[0]-1)] = var[Idx(size[0]-2,size[0]-2)];
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