from red_tide_phase import red_tide_phase
import numpy as np
import matplotlib.pyplot as plt

def red_tide_plot(F, X_Coef, *args):
    varargin = args
    nargin = len(args)
    phaseS, phaseC = red_tide_phase(X_Coef)
    f,s = plt.subplots(2, 1)
    s[0].loglog(F, 0.5*(X_Coef[:,0]**2. + X_Coef[:,1]**2.), '.:')
    s[0].set_xlabel('Frequency (cph)')
    s[0].set_ylabel('Variance contribution (units of time series squared)')
    if nargin==3:
        s[0].set_xscale(varargin[0])
    s[1].semilogx(F, phaseS, '.:', label='\phi_s_i_n')
    s[1].semilogx(F, phaseC, '.:', label='\phi_c_o_s')
    s[1].legend()
    s[1].set_xlabel('Frequency (cph)')
    s[1].set_ylabel('Phase w.r.t. $t=0$')
    s[1].set_ylim(-np.pi, np.pi)
    if nargin==3:
        s[1].set_xscale(varargin[0])
    return