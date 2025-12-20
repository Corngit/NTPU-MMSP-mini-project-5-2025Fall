#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//Complex structure, used for FFT calculation
typedef struct {
    double re, im;
} cpx;

//FFT complex arithmetic
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
static uint32_t read_u32_le(FILE *fp){
    uint8_t b[4];
    fread(b, 1, 4, fp);
    return (uint32_t)b[0]
         | ((uint32_t)b[1] << 8)
         | ((uint32_t)b[2] << 16)
         | ((uint32_t)b[3] << 24);
}

//Read little-endian 16-bit integers
static uint16_t read_u16_le(FILE *fp){
    uint8_t b[2];
    fread(b, 1, 2, fp);
    return (uint16_t)b[0] | ((uint16_t)b[1] << 8);
}

//WAV data structure
typedef struct {
    int sample_rate;      
    int num_channels;     
    int bits_per_sample;  //bits per sample
    int num_samples;      //sample count per channel 
    int16_t *pcm;         //PCM data (interleaved for multi-channel)
} wav_t;


//Read PCM 16-bit mono WAV
static int load_wav_pcm16_mono(const char *path, wav_t *out){
    memset(out, 0, sizeof(*out));
    FILE *fp = fopen(path, "rb");
    if(!fp) return 0;

    char id[4];
    fread(id, 1, 4, fp);              //RIFF
    read_u32_le(fp);
    fread(id, 1, 4, fp);              //WAVE

    uint16_t audio_format = 0;
    uint32_t data_size = 0;
    long data_pos = 0;

    /*Parsing WAV chunk*/
    while(!feof(fp)){
        fread(id, 1, 4, fp);
        uint32_t size = read_u32_le(fp);

        if(memcmp(id, "fmt ", 4) == 0){
            audio_format = read_u16_le(fp);
            out->num_channels = read_u16_le(fp);
            out->sample_rate = read_u32_le(fp);
            read_u32_le(fp);
            read_u16_le(fp);
            out->bits_per_sample = read_u16_le(fp);
            fseek(fp, size - 16, SEEK_CUR);
        }
        else if(memcmp(id, "data", 4) == 0){
            data_pos = ftell(fp);
            data_size = size;
            fseek(fp, size, SEEK_CUR);
        }
        else{
            fseek(fp, size, SEEK_CUR);
        }
    }
    if(audio_format != 1 || out->bits_per_sample != 16){
        fclose(fp);
        return 0;
    }
    /*Read PCM data*/
    fseek(fp, data_pos, SEEK_SET);
    int total_samples = data_size / 2 / out->num_channels;
    out->num_samples = total_samples;
    out->pcm = (int16_t*)malloc(sizeof(int16_t) * total_samples);

    for(int i = 0; i < total_samples; i++){
        fread(&out->pcm[i], sizeof(int16_t), 1, fp);
        for(int ch = 1; ch < out->num_channels; ch++)
            fseek(fp, sizeof(int16_t), SEEK_CUR);
    }
    fclose(fp);
    return 1;
}

//ms to sample
static int ms_to_samples(double ms, int fs){
    return (int)llround(ms * fs / 1000.0);
}


//The specified function type output generates a window
static void build_window(double *w, int P, const char *type){
    if(strcmp(type, "rectangular") == 0){
        for(int n = 0; n < P; n++) w[n] = 1.0;
    }
    else if(strcmp(type, "hamming") == 0){
        for(int n = 0; n < P; n++)
            w[n] = 0.54 - 0.46 * cos(2.0 * M_PI * n / (P - 1));
    }
}


//main(core),STFT and spectrogram calculations
int main(int argc, char** argv){
    if(argc != 7){
        printf("Usage: spectrogram w_size w_type dft_size f_itv wav_in spec_out\n");
        return 1;
    }

    double w_size_ms = atof(argv[1]);
    const char *w_type = argv[2];
    double dft_ms = atof(argv[3]);
    double hop_ms = atof(argv[4]);

    wav_t wav;
    if(!load_wav_pcm16_mono(argv[5], &wav)){
        printf("Failed to load wav file\n");
        return 1;
    }

    int P = ms_to_samples(w_size_ms, wav.sample_rate);
    int N = ms_to_samples(dft_ms, wav.sample_rate);
    int M = ms_to_samples(hop_ms, wav.sample_rate);

    if(P > N || !is_power_of_two(N)){
        printf("Invalid parameters\n");
        return 1;
    }
    double *window = (double*)malloc(sizeof(double) * P);
    cpx *buf = (cpx*)malloc(sizeof(cpx) * N);
    build_window(window, P, w_type);

    FILE *fo = fopen(argv[6], "w");

    for(int s = 0; s < wav.num_samples; s += M){
        for(int n = 0; n < N; n++){
            double x = 0.0;
            if(n < P && s + n < wav.num_samples)
                x = (wav.pcm[s + n] /32768.0) * window[n];
            buf[n].re = x;
            buf[n].im = 0.0;
        }

        fft(buf, N);
    
    //Output single-sided spectrum
    const double eps = 1e-12;

    for(int k = 0; k <= N / 2; k++){
        double mag = sqrt(buf[k].re * buf[k].re + buf[k].im * buf[k].im);

    // normalize by FFT length
        mag /= (double)N;

    // single-sided amplitude correction (except DC and Nyquist)
        if(k != 0 && k != N/2) mag *= 2.0;

        double db = 20.0 * log10(mag + eps);
        fprintf(fo, "%.15f ", db);
   }
    fprintf(fo, "\n");
}


    fclose(fo);
    free(window);
    free(buf);
    free(wav.pcm);
    return 0;

}