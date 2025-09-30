import sigmf
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import scipy.fft as fft

handle = sigmf.sigmffile.fromfile("res/fm_radio_20250920_6msps.sigmf")
handle.read_samples()

plt.specgram(handle[0:10000000], NFFT=8192, Fs=handle.get_global_field(sigmf.SigMFFile.SAMPLE_RATE_KEY))
#plt.plot(np.abs(handle[0:1000]))
#plt.plot(np.abs(fft.fft(handle[0:8192])))
plt.show()

