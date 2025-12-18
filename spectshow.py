import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile

if len(sys.argv) != 4:
    print("Usage: python3 spectshow.py in.wav in.txt out.pdf")
    sys.exit(1)

in_wav = sys.argv[1]
in_txt = sys.argv[2]
out_pdf = sys.argv[3]

## Read wav file
fs, x = wavfile.read(in_wav)
if x.ndim == 2:
    x = x[:, 0]
t = np.arange(len(x)) / fs

## Read spectrogram data
data = np.loadtxt(in_txt)

## Plot spectrogram
fig, ax = plt.subplots(2, 1, figsize=(10, 6))

ax[0].plot(t, x)
ax[0].set_title('Waveform')
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('Amplitude')

## change to dB scale and plot spectrogram
eps = 1e-10
S = data.T      ## here, we assume data is in shape (time_frames, freq_bins) be careful
S_db = 20 * np.log10(np.abs(S) + eps)

vmax = np.max(S_db)
vmin = vmax - 60  

ax[1].imshow(
    S_db,
    aspect='auto',
    origin='lower',
    cmap='gray',     
    vmin=vmin,
    vmax=vmax
)

ax[1].set_title('Spectrogram')
ax[1].set_xlabel('Time Frames')
ax[1].set_ylabel('Frequency Bins')

plt.tight_layout()
plt.savefig(out_pdf)
plt.close(fig)

