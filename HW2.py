import numpy as np
import matplotlib.pyplot as plt

k = int(input("Enter k parameter: "))

N = 2*k + 1
center = k

H = np.zeros(N, dtype=complex)

for m in range(1, center):
    H[m] = -1j

for m in range(center+1, N):
    H[m] = 1j

r_temp = np.fft.ifft(H)
r_temp = np.real(r_temp)

r = np.zeros(N)

for n in range(N):
    index = (n - center) % N
    r[n] = r_temp[index]

n_vals = np.arange(-center, center+1)

NFFT = 2048
w = np.linspace(-np.pi, np.pi, NFFT, endpoint=False)

R = np.zeros(NFFT, dtype=complex)

for i in range(len(n_vals)):
    ni = n_vals[i]
    R += r[i] * np.exp(-1j * w * ni)

freq = w / np.pi

fig, ax = plt.subplots(2,1, figsize=(10,8))

plt.suptitle("Discrete Hilbert Transform Filter")

imagR = np.imag(R)

ax[0].plot(freq, imagR, 'b-', linewidth=2, label='Imag part')
ax[0].plot([-1,0],[1,1],'r--', label='Ideal')
ax[0].plot([0,1],[-1,-1],'r--')

ax[0].axhline(0,color='k',linewidth=0.5)
ax[0].axvline(0,color='k',linestyle=':')

ax[0].set_xlim(-1,1)
ax[0].set_ylim(-1.5,1.5)

ax[0].set_xlabel("Normalized frequency")
ax[0].set_ylabel("Imag part")
ax[0].set_title("Frequency Response")
ax[0].grid(True)
ax[0].legend()

ax[1].stem(n_vals, r, linefmt='b-', markerfmt='bo', basefmt='k-')

ax[1].set_xlabel("n")
ax[1].set_ylabel("r[n]")
ax[1].set_title("Impulse Response")
ax[1].grid(True)

plt.tight_layout()

fname = "hilbert_k" + str(k) + ".png"
plt.savefig(fname)

plt.show()