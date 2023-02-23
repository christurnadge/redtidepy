import numpy as np

def H_make(T, F):
    #if np.shape(T)[1]>np.shape(T)[0]:
    #    T = T.T
    H = np.zeros([len(T), 2*len(F)])
    for i in range(len(F)):
        H[:, 2*i  ] = np.sin(2.*np.pi*F[i]*T)
        H[:, 2*i+1] = np.cos(2.*np.pi*F[i]*T)
    return H