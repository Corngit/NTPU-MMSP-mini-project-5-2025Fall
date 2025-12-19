#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define M_PI 3.14159265358979323846

//複數結構，FFT計算用
typedef struct {
    double re, im;
} cpx;

//FFT計算工具
static cpx c_add(cpx a, cpx b){
    
}

static cpx c_sub(cpx a, cpx b){

}

static cpx c_mul(cpx a, cpx b){

}

//A radix-2 FFT requires an FFT length of 2^k.
static int is_power_of_two(int n){

}

//forward transform
static void fft(cpx *a, int n){
    /*Bit-reversal permutation*/

    /*Butterfly Calculation*/
}

//Reading little-endian 32-bit integers


//Read little-endian 16-bit integers


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



//指定的函數類型輸出產生window



//main(core),STFT and spectrogram calculations
int main(int argc, char** argv){


}