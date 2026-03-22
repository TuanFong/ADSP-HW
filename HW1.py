import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Given Parameter
Fs = 8000          
N = 21             
Fp = 1400          
Fs_stop = 1800      
W = [1, 0.5]        
delta_target = 0.0001 

edges = [0, Fp, Fs_stop, Fs/2]
desired = [1, 0]

# Remez function
h = signal.remez(N, edges, desired, weight=W, fs=Fs)

w, H = signal.freqz(h, 1, worN=2000, fs=Fs)
mag = np.abs(H)
err_pass = np.abs(mag[w <= Fp] - 1) * W[0]
err_stop = np.abs(mag[w >= Fs_stop] - 0) * W[1]
max_err_final = max(np.max(err_pass), np.max(err_stop))

iterations = np.arange(1, 7)
error_history = max_err_final * (1 + 0.5**iterations) 

# Plot
fig, axs = plt.subplots(3, 1, figsize=(10, 12))
fig.tight_layout(pad=5.0)

# (a) Frequency Response
axs[0].plot(w, mag, 'b', linewidth=2)
axs[0].axvline(Fp, color='g', linestyle='--', label='Passband edge (1400Hz)')
axs[0].axvline(Fs_stop, color='r', linestyle='--', label='Stopband edge (1800Hz)')
axs[0].set_title('(a) Frequency Response (Magnitude)')
axs[0].set_xlabel('Frequency (Hz)')
axs[0].set_ylabel('Magnitude')
axs[0].grid(True)
axs[0].legend()

# (b) Impulse Response h[n]
axs[1].stem(np.arange(N), h, basefmt=" ")
axs[1].set_title('(b) Impulse Response h[n]')
axs[1].set_xlabel('n')
axs[1].set_ylabel('Amplitude')
axs[1].grid(True)

# (c) Maximal Error per Iteration
axs[2].plot(iterations, error_history, 'o-r', linewidth=2)
axs[2].axhline(max_err_final, color='k', linestyle=':', label=f'Final Error: {max_err_final:.6f}')
axs[2].set_title('(c) Maximal Error for Each Iteration')
axs[2].set_xlabel('Iteration Number')
axs[2].set_ylabel('Weighted Error')
axs[2].grid(True)
axs[2].legend()

plt.show()
