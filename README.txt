This is the WhoSGlAd method (Farnir et al. 2019) implemented in python.
The WhoSGlAd module can be run as is from tthe terminal or imported to
be used by an external script. When executing from the terminal, the
frequency file to analyse must be supplied as an argument: 
./WhoSGlAd.py <freqfile>

The WhoSGlAd.py script automatically adjusts both the helium and base of
the convection glitches, simultaneously with the smooth contribution to 
theoscillation spectrum. Two plots are produced:
- <prefix>-echelle.pdf: The echelle diagram (nu vs nu%Delta nu) for both
  fitted and reference data. <prefix> being the prefix of the input 
  frequency file.
- <prefix>-glitch.pdf: fThe extracted glitch from the reference data and
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
of the reference frequencies.
