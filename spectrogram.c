#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define M_PI 3.14159265358979323846

//Complex structure, used for FFT calculation
typedef struct {
    double re, im;
} cpx;

//FFT計算工具
static cpx c_add(cpx a, cpx b){
    cpx r = { a.re + b.re, a.im + b.im };
    return r;
}

static cpx c_sub(cpx a, cpx b){
    cpx r = { a.re - b.re, a.im - b.im };
    return r;
}

static cpx c_mul(cpx a, cpx b){
    cpx r = {
        a.re * b.re - a.im * b.im,
        a.re * b.im + a.im * b.re
    };
    return r;
}

//A radix-2 FFT requires an FFT length of 2^k.
static int is_power_of_two(int n){
    return (n > 0) && ((n & (n - 1)) == 0);
}

//forward transform
static void fft(cpx *a, int n){
    /*Bit-reversal permutation*/  
    //a:Input/output complex array
    //n:FFT length (must be a power of 2)
    int j = 0;
    for(int i = 1; i < n; i++){
        int bit = n >> 1;
        while(j & bit){
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if(i < j){
            cpx tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
    }                              
    /*Butterfly Calculation*/
    for(int len = 2; len <= n; len <<= 1){
        double ang = -2.0 * M_PI / (double)len;
        cpx wlen = { cos(ang), sin(ang) };

        for(int i = 0; i < n; i += len){
            cpx w = {1.0, 0.0};
            for(int k = 0; k < len / 2; k++){
                cpx u = a[i + k];
                cpx v = c_mul(a[i + k + len / 2], w);
                a[i + k] = c_add(u, v);
                a[i + k + len / 2] = c_sub(u, v);
                w = c_mul(w, wlen);
            }
        }
    }
}

//Reading little-endian 32-bit integers
static uint32_t read_u32_le(FILE *fp){}

//Read little-endian 16-bit integers
static uint16_t read_u16_le(FILE *fp){}

//WAV data structure
typedef struct {
    int sample_rate;      
    int num_channels;     
    int bits_per_sample;  //每個sample位元數
    int num_samples;      //每個聲道的sample數
    int16_t *pcm;         //單聲道PCM資料
} wav_t;

//Read PCM 16-bit mono WAV
static int load_wav_pcm16_mono(const char *path, wav_t *out){
    /*Parsing WAV chunk*/

    /*Read PCM data*/
}


//毫秒轉sample
static int ms_to_samples(double ms, int fs){}


//指定的函數類型輸出產生window
static void build_window(double *w, int P, const char *type){}


//main(core),STFT and spectrogram calculations
int main(int argc, char** argv){


}