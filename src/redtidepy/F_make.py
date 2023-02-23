import numpy as np

def F_make(*args):
    varargin = args
    nargin = len(args)
    f_m2 =        1./(360./28.9841042 ) # M2
    f_s2 =        1./(360./30.        ) # S2
    f_n2 =        1./(360./28.4397295 ) # N2
    f_k1 =        1./(360./15.0410686 ) # K1
    f_m4 =        1./(360./57.9682084 ) # M4
    f_o1 =        1./(360./13.9430356 ) # O1
    f_m6 =        1./(360./86.9523127 ) # M6
    f_mk3 =       1./(360./44.0251729 ) # MK3
    f_s4 =        1./(360./60.        ) # S4
    f_mn4 =       1./(360./57.4238337 ) # MN4
    f_nu2 =       1./(360./28.5125831 ) # NU2
    f_s6 =        1./(360./90.        ) # S6
    f_mu2 =       1./(360./27.9682084 ) # MU2
    f_2n2 =       1./(360./27.8953548 ) # 2N2
    f_oo1 =       1./(360./16.1391017 ) # OO1
    f_lam2 =      1./(360./29.4556253 ) # LAM2
    f_s1 =        1./(360./15.        ) # S1
    f_m1 =        1./(360./14.4920521 ) # M1
    f_j1 =        1./(360./15.5854433 ) # J1
    f_mm =        1./(360./0.5443747  ) # MM
    f_ssa =       1./(360./0.0821373  ) # SSA
    f_sa =        1./(360./0.0410686  ) # SA
    f_msf =       1./(360./1.0158958  ) # MSF
    f_mf =        1./(360./1.0980331  ) # MF
    f_rho =       1./(360./13.4715145 ) # RHO
    f_q1 =        1./(360./13.3986609 ) # Q1
    f_t2 =        1./(360./29.9589333 ) # T2
    f_r2 =        1./(360./30.0410667 ) # R2
    f_2q1 =       1./(360./12.8542862 ) # 2Q1
    f_p1 =        1./(360./14.9589314 ) # P1
    f_2sm2 =      1./(360./31.0158958 ) # 2SM2
    f_m3 =        1./(360./43.4761563 ) # M3
    f_l2 =        1./(360./29.5284789 ) # L2
    f_2mk3 =      1./(360./42.9271398 ) # 2MK3
    f_k2 =        1./(360./30.0821373 ) # K2
    f_m8 =        1./(360./115.9364169) # M8
    f_ms4 =       1./(360./58.9841042 ) # MS4
    if nargin==9:
        df = varargin[0]
        n_lowNO = varargin[1]
        df_NO = varargin[2]
        n_lowO = varargin[3]
        tide_f = varargin[4]
        sideband_centers = varargin[5]
        n_sidebands = varargin[6]
        df_sidebands = varargin[7]
        inertial = varargin[8]
        if np.shape(n_sidebands)!=np.shape(sideband_centers):
            print('size(n_sidebands) must == size(sideband_centers)')
        Tide_Cell = ['M2','S2','N2','K1','M4','O1','M6','MK3','S4','MN4',
                     'Nu2','S6','MU2','2N2','OO1','Lam2','S1','M1','J1','Mm',
                     'Ssa','Sa','Msf','Mf','Rho','Q1','T2','R2','2Q1','P1',
                     '2SM2','M3','L2','2MK3','K2','M8','MS4']
        Full_Tide_Vec = [f_m2, f_s2, f_n2, f_k1, f_m4, f_o1, f_m6, f_mk3,
                         f_s4, f_mn4, f_nu2, f_s6, f_mu2, f_2n2, f_oo1, 
                         f_lam2, f_s1, f_m1, f_j1, f_mm, f_ssa, f_sa, f_msf,
                         f_mf, f_rho, f_q1, f_t2, f_r2, f_2q1, f_p1, f_2sm2, 
                         f_m3, f_l2, f_2mk3, f_k2, f_m8, f_ms4]
        F_user_defined_tide = []
        F_user_defined_center = []
        if 'tide_f' in globals():
            F_tidal_boolean = np.zeros(np.shape(Tide_Cell))
            F_cusp_boolean = np.zeros(np.shape(Tide_Cell))
            for i in range(len(tide_f)):
                for j in range(len(F_tidal_boolean)):
                    if type(tide_f[i])==str:
                        F_tidal_boolean[j] = eval('Tide_Cell[j]==tide_f[i]')
                    else:
                        F_user_defined_tide = [F_user_defined_tide, tide_f[i]]
            for i in range(len(sideband_centers)):
                for k in range(len(F_cusp_boolean)):
                    if np.isempty(sideband_centers):
                        F_user_defined_center = []
                    elif type(sideband_centers[i]==str):
                        F_cusp_boolean[k] = eval('Tide_Cell[k]==sideband_centers[i]')
                    else:
                        F_user_defined_center = [F_user_defined_center, 
                                                 sideband_centers[i]]
            F_tidal = Full_Tide_Vec(logical(F_tidal_boolean))                  # <<< ??? (logical)
            F_tidal = [F_tidal, F_user_defined_tide]
            F_tidal = np.unique(F_tidal)
            F_cuspcenters = Full_Tide_Vec(logical(F_cusp_boolean))
            F_cuspcenters = [F_cuspcenters, F_user_defined_center]
            F_cuspcenters = np.unique(F_cuspcenters)
        else:
            F_tidal = []
            F_cuspcenters = []
        F_lowNO = range(df_NO, n_lowNO*df_NO, df_NO)
        if n_lowNO==0:
            F_lowO = range(df, n_lowO*df, df)
        else:
            F_lowO = range(F_lowNO[end]+df, F_lowNO[end]+n_lowO*df, df)
        F_cusps = []
        for i in range(len(F_cuspcenters)):
            F_cusps = [F_cusps, 
                       range(F_cuspcenters[i]-n_sidebands[i]*df_sidebands, F_cuspcenters[i]-df_sidebands, df_sidebands),
                       range(F_cuspcenters[i]+df_sidebands, F_cuspcenters[i]+n_sidebands[i]*df_sidebands, df_sidebands)]
        if np.isempty(inertial):
            F_inertial = []
        elif len(inertial)==3:
            f_cor = 2.*1./23.9344699*np.sind(inertial[1])
            df_inertial = inertial[2]
            F_inertial = range(f_cor, f_cor+inertial[3]*df_inertial, df_inertial)
        else:
            print('The input variable "inertial" is formatted incorrectly.')
        F = np.atleast_2d([F_lowNO, F_lowO, F_inertial, F_cusps, F_cuspcenters, F_tidal]).T
        F = np.unique(F)
        if df_sidebands==0 or ~np.isfinite(df_sidebands):
            df_sidebands = df
        diff_F = np.diff(F)
        if np.isempty(F_lowNO):
            if np.sum(diff_F)<0.99*df_sidebands:
                print('In your F vector, you have a frequency difference '+
                      'less than 1/('+str((1./df_sidebands))+' hr) = '+
                      'df_sidebands (or the fundamental frequency '+
                      ' df = 1/(record length) if df_sidebands was given as '+
                      '0, i.e. skipped). This could lead to very high '+
                      'uncertainty at this and nearby frequencies.'+
                      ' The corresponding F indices are:')
                F_Indices = range(len(F))
                print(F_Indices[diff_F<0.99*df_sidebands])
        else:
            if np.sum(diff_F)<0.99*df_NO:
                print('In your F vector, you have a frequency difference '+
                      'less than the lowest you prescribed in df_NO = 1/('+
                      str((1./df_NO))+' hr). This could lead to very high '+
                      'uncertainty at this and nearby frequencies. The '+
                      'corresponding F indices are:')
                F_Indices = range(len(F))
                print(F_Indices[diff_F<0.99*df_NO])
    elif nargin==1:
        NOAA_HARM_CONST = [[1,  f_m2  ], # M2
                           [2,  f_s2  ], # S2
                           [3,  f_n2  ], # N2
                           [4,  f_k1  ], # K1
                           [5,  f_m4  ], # M4
                           [6,  f_o1  ], # O1
                           [7,  f_m6  ], # M6
                           [8,  f_mk3 ], # MK3
                           [9,  f_s4  ], # S4
                           [10, f_mn4 ], # MN4
                           [11, f_nu2 ], # NU2
                           [12, f_s6  ], # S6
                           [13, f_mu2 ], # MU2
                           [14, f_2n2 ], # 2N2
                           [15, f_oo1 ], # OO1
                           [16, f_lam2], # LAM2
                           [17, f_s1  ], # S1
                           [18, f_m1  ], # M1
                           [19, f_j1  ], # J1
                           [20, f_mm  ], # MM
                           [21, f_ssa ], # SSA
                           [22, f_sa  ], # SA
                           [23, f_msf ], # MSF
                           [24, f_mf  ], # MF
                           [25, f_rho ], # RHO
                           [26, f_q1  ], # Q1
                           [27, f_t2  ], # T2
                           [28, f_r2  ], # R2
                           [29, f_2q1 ], # 2Q1
                           [30, f_p1  ], # P1
                           [31, f_2sm2], # 2SM2
                           [32, f_m3  ], # M3
                           [33, f_l2  ], # L2
                           [34, f_2mk3], # 2MK3
                           [35, f_k2  ], # K2
                           [36, f_m8  ], # M8
                           [37, f_ms4 ]] # MS4
        F = np.empty([varargin])*np.nan
        for i in range(varargin):
            F[i] = NOAA_HARM_CONST[NOAA_HARM_CONST[:][1]==i][1]
        F = np.unique(F)
        return F
    else:
        print('Wrong number of inputs.')
