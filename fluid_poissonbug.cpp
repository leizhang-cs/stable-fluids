#include "fluid.h"
#include "pixel.h"


template class Fluid<float,2>;
using pixel = unsigned int;

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
            int val = std::min(U0[i][j].magnitude(), (T)255);
            image_color[i*size[1]+j] = make_pixel(val, val, val);
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
    int i, j;
    XtoIdx(X, i, j);
    U1[i][j] = U0[i][j] + dm;
    std::cout<<U1[i][j]<<std::endl;
}

template<class T, int n>
void Fluid<T,n>::Advect(){
    // U1, U0, dt
    // TraceParticle: method of charactristic
    for(int i=0; i<N[0]; i++){
        for(int j=0; j<N[1]; j++){
            vec X1, X0;
            vec index1(i, j);
            X1 = O + index1*D;
            TraceParticle(X1, X0);
            // interpolate
            int i0, j0;
            XtoIdx(X0, i0, j0);
            U1[i][j] += U0[i0][j0];
        }
    }
    
}

template<class T, int n>
void Fluid<T,n>::Diffuse(){
    // U1, U0, visc, dt 
    // conjugate gradient. FTCS. BTCS???
    // unknown: U1[i][j]
    // TODO: dont consider boundry condition
    T k = visc*dt; // sign -1?
    // U1[0][0] += k*(U0[1][0] - U0[0][0])/(D[0]*D[0]) + (U0[0][1] - U0[0][0])/(D[1]*D[1]);
    // U1[1][0] += k*(U0[2][0] - 2.0f*U0[1][0] + U0[0][0])/(D[0]*D[0]) + (U0[1][1] - U0[1][0])/(D[1]*D[1]);
    // U1[0][1] += k*(U0[1][1] - U0[0][1])/(D[0]*D[0]) + (U0[0][2] - 2.0f*U0[0][1] + U0[0][0])/(D[1]*D[1]);
    for(int i=1; i<N[0]-1; i++){
        for(int j=1; j<N[1]-1; j++){
            // TODO: optimization
            vec ux = (U0[i+1][j] - 2.0f*U0[i][j] + U0[i-1][j])/(D[0]*D[0]);
            vec uy = (U0[i][j+1] - 2.0f*U0[i][j] + U0[i][j-1])/(D[1]*D[1]);
            U1[i][j] += k*(ux + uy);
        }
    }
    boundry_condition(U1);
}

template<class T, int n>
void Fluid<T,n>::Project(){
    // U1, U0, dt
    T dx2 = D[0]*D[0], dy2 = D[1]*D[1], d2 = dx2*dy2;
    
    // divU
    for(int i=1; i<N[0]-1; i++){
        for(int j=1; j<N[1]-1; j++){
            div[i][j] = 0.5f*((U1[i+1][j][0]-U1[i-1][j][0])/D[0] + (U1[i][j+1][1]-U1[i][j-1][1])/D[1]);
            P[i][j] = 0.0f;
        }
    }
    boundry_condition(div);

    // calculate P
    int i0 = size[0]/2, j0 = size[1]/2;
    T prev = 0.0f;
    for(int iter=0; iter<10; iter++){
        prev = P[i0][j0];
        for(int i=1; i<N[0]-1; i++){
            for(int j=1; j<N[0]-1; j++){
                T A = (P[i+1][j]+P[i-1][j])*dy2;
                T B = (P[i][j+1]+P[i][j-1])*dx2;
                P[i][j] = (A+B-div[i][j]*d2)/(2.0f*(dx2+dy2));
            }
        }
        boundry_condition(P);
        std::cout<<P[i0][j0]<<" "<<P[i0][j0]-prev<<std::endl;
    }

    // w4 = w3 - divP
    T max_val = 0;
    for(int i=1; i<N[0]-1; i++){
        for(int j=1; j<N[1]-1; j++){
            vec gradP(0.5f*(P[i+1][j]-P[i-1][j])/D[0], 0.5f*(P[i][j+1]-P[i][j-1])/D[1]);
            U1[i][j] -= gradP;
            max_val = std::max(max_val, U1[i][j].magnitude());
        }
    }
    std::cout<<"U after P:"<<max_val<<std::endl;
    boundry_condition(U1);
}

template<class T, int n>
void Fluid<T,n>::boundry_condition(std::vector<std::vector<T>>& var){
    var[0][0] = var[1][1];
    var[0][1] = var[1][1];
    var[1][0] = var[1][1];
    var[N[0]-1][N[0]-1] = var[N[0]-2][N[0]-2];
    var[N[0]-1][N[0]-2] = var[N[0]-2][N[0]-2];
    var[N[0]-2][N[0]-1] = var[N[0]-2][N[0]-2];
}

template<class T, int n>
void Fluid<T,n>::boundry_condition(std::vector<std::vector<vec>>& var){
    var[0][0] = var[1][1];
    var[0][1] = var[1][1];
    var[1][0] = var[1][1];
    var[N[0]-1][N[0]-1] = var[N[0]-2][N[0]-2];
    var[N[0]-1][N[0]-2] = var[N[0]-2][N[0]-2];
    var[N[0]-2][N[0]-1] = var[N[0]-2][N[0]-2];
}

template<class T, int n>
void Fluid<T,n>::TraceParticle(vec& X1, vec& X0)
{
    // TODO RK
    int i, j;
    XtoIdx(X1, i, j);
    X0 = X1 - dt*U0[i][j];
}

template<class T, int n>
Vec<T,n> Fluid<T,n>::Interpolate(std::vector<std::vector<vec>>& U, vec& X)
{
    // X U
    int i = std::min(X[0]*L[0]/D[0], N[0]-1), j = std::min(X[1]*L[1]/D[1], N[1]-1);
    return U[i][j];
}

template<class T, int n>
inline int Fluid<T,n>::Idx(int i, int j){
    return i*N[1] + j;
}

template<class T, int n>
void Fluid<T,n>::XtoIdx(vec& X, int& i, int& j){
    i = std::min(X[0]*L[0]/D[0], N[0]-1);
    j = std::min(X[1]*L[1]/D[1], N[1]-1);
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