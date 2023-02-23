import numpy as np

def hph(H, P):
    HH = H*H[0,:]
    C = HH*np.diag(P)
    return C