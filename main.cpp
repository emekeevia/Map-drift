#include <iostream>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <cmath>
#include "extra_tools.h"

using namespace std;

void theta(vector<complex<double>> & s_f_1, vector<complex<double>> & s_f_2){

}

typedef complex<double> base;
void fft (vector<complex<double>> & a, bool invert, double scale = 1.0){
    int n = (int) a.size();
    for (int i=1, j=0; i<n; ++i) {
        int bit = n >> 1;
        for (; j>=bit; bit>>=1){
            j -= bit;
        }
        j += bit;
        if (i < j) swap (a[i], a[j]);
    }
    for (int len=2; len<=n; len<<=1) {
        double ang = 2*M_PI*scale/len * (invert ? -1 : 1);
        complex<double> wlen(cos(ang), sin(ang));
        for (int i=0; i<n; i+=len) {
            base w (1);
            for (int j=0; j<len/2; ++j) {
                base u = a[i+j],  v = a[i+j+len/2] * w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }
    if (invert)
        for (int i=0; i<n; ++i) a[i] /= n;
}

void Map_drift(){
    vector<vector<std::complex<double>>> Raw_data = read_file(false, "/home/ivanemekeev/CLionProjects/SAR-data/DIRSIG_input.txt");

    //1)Divide into 2 subapertures
    double k = 0.0, T_a = 0.0;
    vector<vector<std::complex<double>>> s_one_t(Raw_data.begin(), Raw_data.begin() + Raw_data.size()*0.5);
    vector<vector<std::complex<double>>> s_two_t(Raw_data.begin() + Raw_data.size()*0.5, Raw_data.end());
    double gamma_1 = 1 - k * (T_a/2);
    double gamma_2 = 1 + k * (T_a/2);

    //2) Scaled FFT
    // потом обернуть всё в цикл

    vector<std::complex<double>> s_one_f = s_one_t[0];
    vector<std::complex<double>> s_two_f = s_two_t[0];

    fft(s_one_f, true, gamma_1);
    fft(s_two_f, true, gamma_2);


}

int main() {
    Map_drift();
    return 0;
}
