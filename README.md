# WhoSGlAd: **Who**le **S**pectrum and **Gl**itches **Ad**justment

```
                                                                ;
                                                    _           ;;
   __________________________________  _  ______ __| |_________ ;';. _
 ||             _                     | |       / _` |          ;  ;;
 ||____________| '__ ________________ | |   _  | (_| |_________ ; _ ;;
 ||            |  _ \                 | |  / \  \__,_|          ;    ;;
 ||____________| | | | ___ __________ |_| / _ \ _______________ ; __ ;;
 || __        _|_| |_|/ _ \____   ____   / ___ \                ;   ;'
 ||_\ \      / /_____| (_)/ ___| / ___| /_/   \_\ _____________ ;  '___
 ||  \ \ /\ / /       \___\___ \| |  _                     ,;;;,;
 ||___\ V  V /____________ ___) | |_| | __________________ ;;;;;; _____
       \_/\_/             /____/ \____|                    `;;;;'
```

This is the WhoSGlAd method (Farnir et al. 2019) implemented in python.
> :memo: Note:
> A fortran implementation exists and can be provided on request.

## Usage
The WhoSGlAd module can be run as is from the terminal or imported to be used by an external script. When executing from the terminal, the frequency file to analyse must be supplied as an argument: 
`./WhoSGlAd.py <freqfile>`

## Description of the files:

#### **WhoSGlAd.py**: Main code
The WhoSGlAd.py script automatically adjusts both the helium and base of the convection zone glitches, simultaneously with the smooth contribution to the oscillation spectrum. 

Two plots are produced:
- \<prefix\>-echelle.pdf: The echelle diagram (nu vs nu%Delta nu) for both   fitted and reference data. \<prefix\> being the prefix of the input frequency file.
- \<prefix\>-glitch.pdf: The extracted glitch from the reference data and   a 'continuous' representation of the fitted glitch functions. <prefix> being the prefix of the input frequency file.

Two output text files are also produced:
- \<prefix\>-fit.txt: The fitted and reference frequencies. \<prefix\> being the prefix of the input frequency file.
- \<prefix\>-indicators.txt: The computed seismic indicators defined with the WhoSGlAd method (Farnir et al. 2019). \<prefix\> being the prefix of the input frequency file.

#### **pltutils.py**: Handy plooting functions and definitions
The script pltutils.py defines some handy functions and the selection of colors to use. It is not essential and may easily be replaced.

#### **WhoSGlAd_cfg.py**: Define WhoSGlAd's behavior
- *T_He and T_CZ*: float: The WhoSGlAd_cfg_py file contains the dimensionless depths of the helium and convection zone glitches needed by WhoSGlAd to carry an adjustment of the reference frequencies. These are defined in Farnir et al. 2019 as   
  T = tau*Delta,   
  where tau is the acoustic depth of the considered glitch, 
  Delta = 1/(2*R_ac) with R_ac the acoustic radius of the star.    
  A good approximation of Delta is the mean large separation (which is also provided by WhoSGlAd). We see that a value for tau is needed but we do not have it in the case of observations. Nevertheless, in the case of the helium glitch it is easily circumvented following Fanir et al. 2023 providing an estimate in terms of observables.
- *T_est*: bool: Whether to use the automated estimated (Farnir et al. 2023) of T_He.
- *idxlist*: dict: WhoSGlAd_cfg_py also specifies the ordering of the input frequency file. It is provided as a dictionnary. Only 4 keywords are recognised:
  - *n*: the radial order.
  - *l*: the spherical degree.
  - *freq*: the frequency value, assumed to be in micro Hz.
  - *sigma*: the frequency uncertainty, assumed to be in micro Hz.
- *target_ln*: [(int,int)]: array of l (spherical degree) and n (radial order) values to be used in the adjustment (crucial to the definition of the basis)
- *plot*: bool: Whether to produce plots
- *save_plots*: bool: Whether to save plots
- *show_plots*: bool: Whether to show plots
- *save_coefs*: bool: Whether to save akl (projection) coefficients
