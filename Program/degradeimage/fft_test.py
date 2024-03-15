from scipy.fft import fft, ifft, fftshift
from numpy import log10
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import windows, find_peaks


# fig, ax = plt.subplots(2, figsize = (12,10))
# rect_100pts = windows.boxcar(50)
# ax[0].plot(rect_100pts)
# ax[0].set_title('Rect window')
# ax[0].set_xlabel('Frequency')
# ax[0].set_ylabel('Magnitude (dB)')

# A = fft(rect_100pts, 128)
# freq = np.linspace(0, 1000, len(A))
# response = 20 * log10(np.abs(fftshift(A)))
# ax[1].plot(freq, response, label = 'lower sampling')
# ax[1].set_title('Frequency response')
# ax[1].set_xlabel('Frequency')
# ax[1].set_ylabel('Magnitude (dB)')

# A = fft(rect_100pts, 1024)
# freq = np.linspace(0, 1000, len(A))
# response = 20 * log10(np.abs(fftshift(A)))
# ax[1].plot(freq, response, label = 'higher sampling')
# ax[1].legend()
# plt.tight_layout()
# # plt.savefig('degradeimage/samplingrates.png')
# plt.show()




# fig, ax = plt.subplots(2, figsize = (12,7))
# hamming_100pts = windows.hamming(100)
# A = fft(hamming_100pts, 2048)
# freq = np.linspace(0, 1000, len(A))
# response = 20 * log10(np.abs(fftshift(A)))
# ax[0].plot(freq, response)
# ax[0].set_title('Hamming with shorter look time')
# ax[0].set_xlabel('Frequency')
# ax[0].set_ylabel('Magnitude (dB)')
# peaks, _ = find_peaks(response)
# firstpeak = peaks[np.argmax(response[peaks])]
# magnitude_3db = response[firstpeak] - 3
# secondpeak = firstpeak - 1
# peakdiff = abs(response[firstpeak] - response[secondpeak])
# print(magnitude_3db)
# cut = np.zeros(len(response))
# cut.fill(magnitude_3db)
# ax[0].plot(cut)
# ax[0].set(xlim = [0, 1000], ylim = [-80, 45])
# idx = np.argwhere(np.diff(np.sign(cut - response))).flatten()
# bandwidth_3db = abs(idx[1] - idx[0])
# ax[0].annotate(f'bandwidth is {bandwidth_3db}Hz', xy = (5,34), size = 8)


# hamming_200pts = windows.hamming(200)
# print(len(hamming_100pts), len(hamming_200pts))
# A = fft(hamming_200pts, 2048)
# freq = np.linspace(0, 1000, len(A))
# response = 20 * log10(np.abs(fftshift(A)))
# ax[1].plot(freq, response)
# ax[1].set_title('Hamming with longer look time')
# ax[1].set_xlabel('Frequency')
# ax[1].set_ylabel('Magnitude (dB)')
# peaks, _ = find_peaks(response)
# firstpeak = peaks[np.argmax(response[peaks])]
# magnitude_3db = response[firstpeak] - 3
# print(magnitude_3db)
# cut = np.zeros(len(response))
# cut.fill(magnitude_3db)
# ax[1].plot(cut)
# ax[1].set(xlim = [0, 1000], ylim = [-80, 45])
# ax[1].annotate('Longer look time results in higher resolution, better ability to resolve 2 distinct objects close together', size = 7, xy = (20,5))
# idx = np.argwhere(np.diff(np.sign(cut - response))).flatten()
# bandwidth_3db = abs(idx[1] - idx[0])
# ax[1].annotate(f'bandwidth is {bandwidth_3db}Hz', xy = (5,30), size = 8)

# plt.tight_layout()
# # plt.savefig('degradeimage/looktime.png')
# plt.show()



fig, ax = plt.subplots(2,2, figsize = (12, 7))

fig.suptitle('Frequency response of Different Windows')
rect_50pts = windows.boxcar(50)
A = fft(rect_50pts, 1024)
freq = np.linspace(0, 1000, len(A))
response_rect = 20 * log10(np.abs(fftshift(A)))
ax[0,0].plot(freq, response_rect)
ax[0,0].set_title('Rectangle')
ax[0,0].set_xlabel('Frequency (Hz)')
ax[0,0].set_ylabel('Magnitude (dB)')
peaks, _ = find_peaks(response_rect)
firstpeak = peaks[np.argmax(response_rect[peaks])]
magnitude_3db = response_rect[firstpeak] - 3
peakdiff = abs(response_rect[firstpeak] - response_rect[peaks[np.argmax(response_rect[peaks])-1]])
cut = np.zeros(len(response_rect))
cut.fill(magnitude_3db)
ax[0,0].plot(cut)
idx = np.argwhere(np.diff(np.sign(cut - response_rect))).flatten()
bandwidth_3db = abs(idx[1] - idx[0])
ax[0,0].annotate(f'bandwidth is {bandwidth_3db}Hz', xy = (0,32), size = 8)
ax[0,0].annotate(f'magnitude diff between first and second\npeak is {peakdiff}dB', xy = (0, 20), size = 7)

taylor_50pts = windows.taylor(50)
A = fft(taylor_50pts, 1024)
freq = np.linspace(0, 1000, len(A))
response_taylor = 20 * log10(np.abs(fftshift(A)))
ax[0,1].plot(freq, response_taylor)
ax[0,1].set_title('Taylor')
ax[0,1].set_xlabel('Frequency (Hz)')
ax[0,1].set_ylabel('Magnitude (dB)')
peaks, _ = find_peaks(response_taylor)
firstpeak = peaks[np.argmax(response_taylor[peaks])]
magnitude_3db = response_taylor[firstpeak] - 3
peakdiff = abs(response_taylor[firstpeak] - response_taylor[peaks[np.argmax(response_taylor[peaks])-1]])
cut = np.zeros(len(response_taylor))
cut.fill(magnitude_3db)
ax[0,1].plot(cut)
idx = np.argwhere(np.diff(np.sign(cut - response_taylor))).flatten()
bandwidth_3db = abs(idx[1] - idx[0])
# ax[0,1].plot(idx, response_taylor[idx], 'x')
ax[0,1].annotate(f'bandwidth is {bandwidth_3db}Hz', xy = (0,29), size = 8)
ax[0,1].annotate(f'magnitude diff between first and second\npeak is {peakdiff}dB', xy = (0, 15), size = 7)

hamming_50pts = windows.hamming(50)
A = fft(hamming_50pts, 1024)
freq = np.linspace(0, 1000, len(A))
response_hamming = 20 * log10(np.abs(fftshift(A)))
ax[1,0].plot(freq, response_hamming)
ax[1,0].set_title('Hamming')
ax[1,0].set_xlabel('Frequency (Hz)')
ax[1,0].set_ylabel('Magnitude (dB)')
peaks, _ = find_peaks(response_hamming)
firstpeak = peaks[np.argmax(response_hamming[peaks])]
magnitude_3db = response_hamming[firstpeak] - 3
peakdiff = abs(response_hamming[firstpeak] - response_hamming[peaks[np.argmax(response_hamming[peaks])-1]])
cut = np.zeros(len(response_hamming))
cut.fill(magnitude_3db)
ax[1,0].plot(cut)
idx = np.argwhere(np.diff(np.sign(cut - response_hamming))).flatten()
bandwidth_3db = abs(idx[1] - idx[0])
ax[1,0].annotate(f'bandwidth is {bandwidth_3db}Hz', xy = (0,28), size = 8)
ax[1,0].annotate(f'magnitude diff between first and second\npeak is {peakdiff}dB', xy = (0, 13), size = 7)

hann_50pts = windows.hann(50)
A = fft(hann_50pts, 1024)
freq = np.linspace(0, 1000, len(A))
response_hann = 20 * log10(np.abs(fftshift(A)))
ax[1,1].plot(freq, response_hann)
ax[1,1].set_title('Hann')
ax[1,1].set_xlabel('Frequency (Hz)')
ax[1,1].set_ylabel('Magnitude (dB)')
peaks, _ = find_peaks(response_hann)
firstpeak = peaks[np.argmax(response_hann[peaks])]
magnitude_3db = response_hann[firstpeak] - 3
peakdiff = abs(response_hann[firstpeak] - response_hann[peaks[np.argmax(response_hann[peaks])-1]])
cut = np.zeros(len(response_hann))
cut.fill(magnitude_3db)
ax[1,1].plot(cut)
idx = np.argwhere(np.diff(np.sign(cut - response_hann))).flatten()
bandwidth_3db = abs(idx[1] - idx[0])
ax[1,1].annotate(f'bandwidth is {bandwidth_3db}Hz', xy = (0,28), size = 8)
ax[1,1].annotate(f'magnitude diff between first and second\npeak is {peakdiff}dB', xy = (0, 10), size = 7)

plt.tight_layout()
# plt.savefig('degradeimage/freqresponsecomparison.png')
plt.show()