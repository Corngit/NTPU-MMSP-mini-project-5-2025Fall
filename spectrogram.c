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



//WAV檔讀取
typedef struct {
    int sample_rate;      
    int num_channels;     
    int bits_per_sample;  //每個sample位元數
    int num_samples;      //每個聲道的sample數
    int16_t *pcm;         //單聲道PCM資料
} wav_t;


//毫秒轉sample



//指定的函數類型輸出產生window



//main(core)
int main(int argc, char** argv){}