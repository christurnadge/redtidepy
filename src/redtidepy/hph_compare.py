import numpy as np
import matplotlib.pyplot as plt
from hph import hph

def hph_compare(H, P, X):
    HPH = hph(H, P)
    xcovX = np.cov(X, bias=True)
    xcovX = np.flip(np.fft.fftshift(xcovX))
    xcovX = xcovX[1:int((len(xcovX)+1)/2)]
    f,s = plt.subplots()
    s.plot(HPH, '.-', label='HPH$^T$')
    s.plot(xcovX, '.-', label='covariance')
    s.legend()
    s.set_xlabel('lag')
    return