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
static int load_wav_pcm16_mono(const char *path, wav_t *out) {
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open wav: %s\n", path);
        return -1;
    }

    // RIFF header
    char id[5] = {0};
    if (fread(id, 1, 4, fp) != 4 || (id[0]!='R'||id[1]!='I'||id[2]!='F'||id[3]!='F')) {
        fprintf(stderr, "ERROR: not RIFF\n");
        fclose(fp);
        return -1;
    }
    (void)read_u32_le(fp); // chunk size
    if (fread(id, 1, 4, fp) != 4 || (id[0]!='W'||id[1]!='A'||id[2]!='V'||id[3]!='E')) {
        fprintf(stderr, "ERROR: not WAVE\n");
        fclose(fp);
        return -1;
    }

    // Find "fmt " and "data"
    uint16_t audioFormat = 0, numChannels = 0, bitsPerSample = 0;
    uint32_t sampleRate = 0;
    uint32_t dataSize = 0;
    long dataPos = -1;

    while (!feof(fp)) {
        if (fread(id, 1, 4, fp) != 4) break;
        uint32_t chunkSize = read_u32_le(fp);

        if (id[0]=='f' && id[1]=='m' && id[2]=='t' && id[3]==' ') {
            audioFormat  = read_u16_le(fp);
            numChannels  = read_u16_le(fp);
            sampleRate   = read_u32_le(fp);
            (void)read_u32_le(fp); // byteRate
            (void)read_u16_le(fp); // blockAlign
            bitsPerSample = read_u16_le(fp);

            // skip rest of fmt if any
            uint32_t remain = (chunkSize > 16) ? (chunkSize - 16) : 0;
            if (remain) fseek(fp, (long)remain, SEEK_CUR);

        } else if (id[0]=='d' && id[1]=='a' && id[2]=='t' && id[3]=='a') {
            dataPos = ftell(fp);
            dataSize = chunkSize;
            fseek(fp, (long)chunkSize, SEEK_CUR);

        } else {
            // skip unknown chunk (pad to even)
            fseek(fp, (long)chunkSize, SEEK_CUR);
        }

        if (chunkSize & 1) fseek(fp, 1, SEEK_CUR); // padding
    }

    if (audioFormat != 1) {
        fprintf(stderr, "ERROR: only PCM supported (audioFormat=%u)\n", audioFormat);
        fclose(fp);
        return -1;
    }
    if (!(numChannels == 1 || numChannels == 2)) {
        fprintf(stderr, "ERROR: only mono/stereo supported (channels=%u)\n", numChannels);
        fclose(fp);
        return -1;
    }
    if (bitsPerSample != 16) {
        fprintf(stderr, "ERROR: only 16-bit supported (bits=%u)\n", bitsPerSample);
        fclose(fp);
        return -1;
    }
    if (dataPos < 0 || dataSize == 0) {
        fprintf(stderr, "ERROR: missing data chunk\n");
        fclose(fp);
        return -1;
    }

    // number of mono samples after downmix
    const int bytesPerSample = 2;
    int totalFrames = (int)(dataSize / (numChannels * bytesPerSample)); // frames = samples per channel
    int16_t *mono = (int16_t*)malloc((size_t)totalFrames * sizeof(int16_t));
    if (!mono) {
        fprintf(stderr, "ERROR: malloc failed\n");
        fclose(fp);
        return -1;
    }

    // read data and downmix
    fseek(fp, dataPos, SEEK_SET);

    for (int i = 0; i < totalFrames; i++) {
        int16_t s0 = 0, s1 = 0;
        if (fread(&s0, sizeof(int16_t), 1, fp) != 1) { free(mono); fclose(fp); return -1; }
        if (numChannels == 2) {
            if (fread(&s1, sizeof(int16_t), 1, fp) != 1) { free(mono); fclose(fp); return -1; }

            // average with 32-bit accumulator to avoid overflow
            int32_t sum = (int32_t)s0 + (int32_t)s1;
            mono[i] = (int16_t)(sum / 2);
        } else {
            mono[i] = s0;
        }
    }

    fclose(fp);

    out->sample_rate = (int)sampleRate;
    out->num_samples = totalFrames;
    out->pcm = mono;
    return 0;
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
            if(s + n < wav.num_samples)
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