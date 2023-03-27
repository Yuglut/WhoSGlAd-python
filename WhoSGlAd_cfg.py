#!/usr/bin/env python3
# Configuration file for the main script WhoSGlAd.py defining a few 
# quantities

# Dimensionless helium glitch acoustic depth
T_He = 0.75591E-01
# Dimensionless base of the convection zone glitch acoustic depth
T_CZ = 0.30240E+00
# Whether to estimate T_He via linear relation (Farnir et al. 2023)
T_Est = True 
# Ordering of the columns of the frequency file. The 4 recognised
# keywords are: 'n' the radial orde, 'l' the sperical degree,
# 'freq' the frequency value (assumed to be in micro Hz), 'sigma'
# the uncertainties on frequencies (also in micro Hz). If the number
# of columns in the input file exceeds that of idxlist, the unlist 
# columns won't be considered by WhoSGlAd
idxlist={'l':0,'n':1,'freq':2,'sigma':3}
# A list of modes (l,n) to use for the adjustment
# Typical 16Cyg A
target_ln=[(0,n) for n in range(13,26)] + [(1,n) for n in range(13,25)]\
          +[(2,n) for n in range(12,25)]+[(3,n) for n in range(15,23)]
plot=True # Whether to produce plots
save_plots=True # Whether to save plots
show_plots=True # Whether to show plots
save_coefs=True # Whether to save akl coefficients
