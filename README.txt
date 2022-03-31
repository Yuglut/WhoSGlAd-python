This is the WhoSGlAd method (Farnir et al. 2019) implemented in python.
The WhoSGlAd module can be run as is from tthe terminal or imported to
be used by an external script. When executing from the terminal, the
frequency file to analyse must be supplied as an argument: 
./WhoSGlAd.py <freqfile>

The WhoSGlAd.py script automatically adjusts both the helium and base of
the convection zone glitches, simultaneously with the smooth 
contribution to the oscillation spectrum. Two plots are produced:
- <prefix>-echelle.pdf: The echelle diagram (nu vs nu%Delta nu) for both
  fitted and reference data. <prefix> being the prefix of the input 
  frequency file.
- <prefix>-glitch.pdf: The extracted glitch from the reference data and
  a 'continuous' representation of the fitted glitch functions. <prefix> 
  being the prefix of the input frequency file.
Two output text files are also produced:
- <prefix>-fit.txt: The fitted and reference frequencies. <prefix> being
  the prefix of the input frequency file.
- <prefix>-indicators.txt: The computed seismic indicators defined with
  the WhoSGlAd method (Farnir et al. 2019). <prefix> being the prefix of
  the input frequency file.

The script pltutils.py defines some handy functions and the selection of
colors to use. It is not essential and may easily be replaced.

The WhoSGlAd_cfg_py file contains the dimensionless depths of the helium
and convection zone glitches needed by WhoSGlAd to carry an adjustment
of the reference frequencies. These are defined in Farnir et al. 2019 as
T = tau*Delta, 
where tau is the acoustic depth of the considered glitch,
Delta = 1/(2*R_ac) with R_ac the acoustic radius of the star. 
A good approximation of Delta is the mean large separation (which is 
also provided by WhoSGlAd). We see that a value for tau is needed but we
do not have it in the case of observations. Nevertheless, in the case of
the helium glitch it is easily circumvented thanks to Angelo Valentino's
master thesis work providing an estimate in terms of observables. (not
yet implemented in the current version)
WhoSGlAd_cfg_py also specifies the ordering of the input frequncy file.
It is provided as a dictionnary. Only 4 keywords are recognised:
- 'n': the radial order.
- 'l': the spherical degree.
- 'freq': the frequency value, assumed to be in micro Hz.
- 'sigma': the frequency uncertainty, assumed to be in micro Hz.
