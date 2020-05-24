#include <iostream>
#include <fftw.h>

int main(){
    int N = 10;
    fftw_complex in[N], out[N];
    fftw_plan p;
    
    p = fftw_create_plan(N, FFTW_FORWARD, FFTW_ESTIMATE);
    
    fftw_one(p, in, out);
    
    fftw_destroy_plan(p);  
}