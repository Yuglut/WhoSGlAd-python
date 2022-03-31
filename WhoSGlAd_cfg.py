#!/usr/bin/env python3
# Configuration file for the main script WhoSGlAd.py defining a few 
# quantities

# Dimensionless helium glitch acoustic depth
T_He = 9.0488329068490900E-002
# Dimensionless base of the convection zone glitch acoustic depth
T_CZ = 0.31753506074216359
# Ordering of the columns of the frequency file. The 4 recognised
# keywords are: 'n' the radial orde, 'l' the sperical degree,
# 'freq' the frequency value (assumed to be in micro Hz), 'sigma'
# the uncertainties on frequencies (also in micro Hz)
idxlist={'l':0,'n':1,'freq':2,'sigma':3}
