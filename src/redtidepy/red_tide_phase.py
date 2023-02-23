import numpy as np

def red_tide_phase(X_Coef, *args):
    varargin = args
    nargin = len(args)+1
    if nargin==1:
        phaseSIN = np.arctan2( X_Coef[:,1], X_Coef[:,0])
        phaseCOS = np.arctan2(-X_Coef[:,0], X_Coef[:,1])
        return phaseSIN, phaseCOS 
    elif nargin==2:
        F = varargin[0,0]
        Dateformat = varargin[0,2]
        T0 = varargin[0,2]
        T0_new = varargin[0,3]
        t0 = datenum(T0, Dateformat)                                           # <<< ??? (datenum)
        t0_new = datenum(T0_new, Dateformat)                                   # <<< ??? (datenum)
        dt_hours = (t0-t0_new)*24.
        phaseSIN = np.mod(np.arctan2( X_Coef[:,1],X_Coef[:,0])+dt_hours*2.*
                                      np.pi*F, 2.*np.pi)
        phaseCOS = np.mod(np.arctan2(-X_Coef[:,0],X_Coef[:,1])+dt_hours*2.*
                                      np.pi*F, 2.*np.pi)
        return phaseSIN, phaseCOS 
    else:
        print('Incorrect number of inputs, please see documentation.')
