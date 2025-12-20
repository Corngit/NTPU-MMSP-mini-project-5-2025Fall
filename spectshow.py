import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile

if len(sys.argv) != 4:
    print("Usage: python spectshow.py input.wav input.txt output.pdf")
    sys.exit(1)

wav_path = sys.argv[1]
txt_path = sys.argv[2]
pdf_path = sys.argv[3]


fs, x = wavfile.read(wav_path)
if x.ndim == 2:
    x = x[:, 0]   # mono

t = np.arange(len(x)) / fs

data = np.loadtxt(txt_path)

if data.size == 0:
    raise RuntimeError("Spectrogram txt is empty")


S_db = data.T   


vmax = np.nanpercentile(S_db, 99)
vmin = vmax - 60


fig, ax = plt.subplots(2, 1, figsize=(10, 6))

# Waveform
ax[0].plot(t, x)
ax[0].set_title("Waveform")
ax[0].set_xlabel("Time [s]")
ax[0].set_ylabel("Amplitude")

# Spectrogram
im = ax[1].imshow(
    S_db,
    aspect="auto",
    origin="lower",
    vmin=vmin,
    vmax=vmax
)
ax[1].set_title("Spectrogram")
ax[1].set_xlabel("Time Frames")
ax[1].set_ylabel("Frequency Bins")

fig.colorbar(im, ax=ax[1], label="Magnitude (dB)")

plt.tight_layout()
plt.savefig(pdf_path)
plt.close()
