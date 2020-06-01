#include <iostream>
#include <fftw3.h>

using namespace std;

// g++ -o fft_test -std=c++14 fft_test.cpp -lfftw3 -lm

int main(){
    int N = 10;
    
    // fftw_complex* in = new fftw_complex[N]();
    // fftw_complex* out = new fftw_complex[N]();
    double* in = new double[N]();
    double* out = new double[N]();
    fftw_plan p;
    
    cout<<"O"<<endl;
    for(int i=0; i<N/2; i++){
        in[i] = rand()/(RAND_MAX+1.0);
        //in[i][1] = rand()/(RAND_MAX+1.0);
    }
    for(int i=0; i<N; i++){
        cout<<in[i]<<endl;
    }

    p = fftw_plan_r2r_1d(N, in, out, FFTW_R2HC, FFTW_ESTIMATE);
    
    fftw_execute(p);

    cout<<"F"<<endl;
    for(int i=0; i<N; i++){
        cout<<out[i]<<endl;
    }
    
    double res[N];

    fftw_plan pb = fftw_plan_r2r_1d(N, out, res, FFTW_HC2R, FFTW_ESTIMATE);

    fftw_execute(pb);

    cout<<"B"<<endl;
    for(int i=0; i<N; i++){
        cout<<res[i]/N<<endl;
    }
    
    fftw_destroy_plan(p);
    fftw_destroy_plan(pb);
    
    delete []in;
    delete []out;

    return 0;
}