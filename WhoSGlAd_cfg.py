#!/usr/bin/env python3
# Configuration file for the main script WhoSGlAd.py defining a few 
# quantities

# Dimensionless helium glitch acoustic depth
T_He = 0.75591E-01 # 16CygA
T_He = 9.6282790077356989E-002 # Gemma
#T_He = 0.7406987696896E-01 # Doris
# Dimensionless base of the convection zone glitch acoustic depth
T_CZ = 0.30240E+00 # 16CygA
T_CZ = 0.22706532878684116 # Gemma
#T_CZ = 0.3207150078947E+00 # Doris
# Whether to estimate T_He via linear relation (Farnir et al. 2023)
T_Est = True 
# Ordering of the columns of the frequency file. The 4 recognised
# keywords are: 'n' the radial order, 'l' the sperical degree,
# 'freq' the frequency value (assumed to be in micro Hz), 'sigma'
# the uncertainties on frequencies (also in micro Hz). If the number
# of columns in the input file exceeds that of idxlist, the unlist 
# columns won't be considered by WhoSGlAd
idxlist={'l':0,'n':1,'freq':2,'sigma':3}
# A list of modes (l,n) to use for the adjustment
# Typical 16Cyg A
#target_ln=[(0,n) for n in range(13,26)] + [(1,n) for n in range(13,25)]\
#          +[(2,n) for n in range(12,25)]+[(3,n) for n in range(15,23)]
# Gemma
target_ln=[(0,n) for n in range(10,23)]
# KIC6679371
#target_ln=[(0,n) for n in range(9,26)] + [(1,n) for n in range(9,26)]\
#          +[(2,n) for n in range(12,24)]
# KIC10162436
#target_ln=[(0,n) for n in range(11,26)] + [(1,n) for n in range(8,26)]\
#          +[(2,n) for n in range(11,24)]
# Ok for SunAsAStar Solar Cycle
#target_ln=[(0,n) for n in range(14,26)] + [(1,n) for n in range(13,25)]\
#          +[(2,n) for n in range(13,25)]
# High frequencies SunAsAStar Solar Cycle
#target_ln=[(0,n) for n in range(19,26)] + [(1,n) for n in range(18,25)]\
#          +[(2,n) for n in range(18,25)]
# Doris, Magnetic Cycle
#target_ln=[(0,n) for n in range(18,27)] + [(1,n) for n in range(17,27)]\
#         +[(2,n) for n in range(18,24)]
plot=True # Whether to produce plots
save_plots=True # Whether to save plots
show_plots=True # Whether to show plots
save_coefs=False # Whether to save akl coefficients
