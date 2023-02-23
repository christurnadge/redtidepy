import numpy as np
import scipy as sp

def P_make(f, S, F, *args):
    varargin = args
    nargin = len(args)
    if nargin==3:
        SamplePeriod = 1.
        L = 1./np.min(np.diff(F))
        INTERP_METHOD = 'loglinear'
    elif nargin==5:
        SamplePeriod = varargin[0]
        L = varargin[1]
        INTERP_METHOD = 'loglinear'
    elif nargin==6:
        SamplePeriod = varargin[0]
        L = varargin[1]
        INTERP_METHOD = varargin[2]
    else:
        print('3 or 5 or 6 inputs expected, given '+str(nargin)+'. See '+
              'documentation.')
    f_Ny = 1./(2.*SamplePeriod)
    df = 1./L
    F_ = F
    if F_[1]>df:
        F_ = [df:F_]
    if F_[-1]<f_Ny:
        F_ = [F_:f_Ny]
    diff_F = np.diff(F_)
    F_unmodeled = []
    for i in range(len(diff_F)):
        if diff_F[i]>2.*df:
            F_unmodeled = [F_unmodeled, F_[i]+df:df:F_[i+1]-df]                # <<< ??? (syntax)
    F_full = [F, F_unmodeled.T]                                                # <<< ??? (concatenate instead?)
    F_full = np.unique(F_full)
    if INTERP_METHOD=='loglinear':
        Sxx_Prior = np.exp(np.interp1([f;(f_Ny+df)], np.log([S;S(end)]),       # <<< ??? (syntax)
                                       F_full, 'linear'))
    else:
        Sxx_Prior = np.interp1([f;(f_Ny+df)], [S;S(end)], F_full,              # <<< ??? (syntax)
                                INTERP_METHOD)
    Sxx_Prior[Sxx_Prior<=0] = np.min(Sxx_Prior[Sxx_Prior>0])
    if np.isinf(S[1]):
        print('Your averaged spectral estimates begin with a non-finite '+
              'value, which absolutely should not be the case.')
    Sxx_Prior[np.isnan(Sxx_Prior)] = S[1]
    diffF_full = np.diff(F_full)
        diffF_full = [diffF_full[1]; diffF_full];                              # <<< ??? (syntax)
    Sxx_Prior = Sxx_Prior*(diffF_full/np.median(np.diff(f)))
    Sxx_ModelPrior = Sxx_Prior(np.in1d(F_full, F))
    P_diag = np.zeros([2*len(Sxx_ModelPrior), 1])
    P_diag[0:1:2*len(Sxx_ModelPrior)] = Sxx_ModelPrior                         # <<< ??? (syntax) 
    P_diag[1:1:2*len(Sxx_ModelPrior)] = Sxx_ModelPrior                         # <<< ??? (syntax)
    P = sp.sparse.spdiags(P_diag, 0, len(P_diag), len(P_diag))
    return P