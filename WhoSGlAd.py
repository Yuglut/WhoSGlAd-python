#!/usr/bin/env python3
import numpy as np
import os
import sys
import math
import WhoSGlAd_cfg as config
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pltutils as pltu

class Mode:
  """
  A class defining oscilation modes
  """
  def __init__(self,_n,_l,_value,_sigma):
    """
    Initialises the Mode object.

    :param _n: radial order.
    :type _n: int

    :param _l: spherical degree.
    :type _n: int

    :param _value: frequency (expected to be in muHz).
    :type _value: float

    :param _sigma: uncertainty on frequency (expected to be in muHz).
    :type _sigma: float
    """
    self.n     = _n
    self.l     = _l
    self.value = _value
    self.sigma = _sigma

  def print_me(self):
    """
    Print mode
    """
    print("l={:d}, n={:d}, freq={:.3f} +- {:.2e} (muHz)".format(int(self.l), int(self.n), self.value, self.sigma))
    
def get_freq(freqfile,idxlist={'l':0,'n':1,'freq':2,'sigma':3}):
  """
  Reads the frequency file <freqfile> and returns a vector of modes of 
  class Mode <modes>.

  :param idxlist: a list containing the ordering of nu, n, l and sig in
                  the reference file
  :type idxlist: int
  """
  if not(os.path.exists(freqfile)):
    sys.exit('{:s} file missing'.format(freqfile))
  with open(freqfile,'r') as f:
    lines = f.readlines()
    modes = [Mode(float(line.split()[idxlist['n']]),float(line.split()[idxlist['l']]),float(line.split()[idxlist['freq']]),float(line.split()[idxlist['sigma']])) for line in lines if line.lstrip()[0]!='#']
  return modes

class Indicator:
  """ 
  Seismic indicator class containing its name, value and uncertainty
  """

  def __init__(self,value=None,sigma=None):
    self.value = value
    self.sigma = sigma

class Seismic:
  """
  Class of seismic constraints
  """
  # Default list of seismic indicators (to enforce their ordering)
  indic_list=['Dnu','Dnu0','Dnu1','Dnu2','Dnu3','r01','r02','r03','Delta01','Delta02','Delta03','eps','eps0','eps1','eps2','eps3','AHe','ACZ']

  def __init__(self,_indic_list=indic_list):
    self.modes      = [] # Set of 'observed' modes
    self.fit        = [] # Set of fitted modes
    self.Rkkinv     = None # Inverse R matrix of tranformation
    self.qkinv      = None # Normalised basis vector
    self.indicators = {key:Indicator() for key in _indic_list} # Computed indicators, stored in a dict.

  def find_l_list(self, l_targets, npowers=2):
      """
      Find a list of l values with the following properties:
      
      - each l value only occurs once
      - each l value given in the parameter ``l_targets`` is
        in the result l list, except if there is less than npowers
        modes with this l value
      - if the parameter ``l_targets`` is ``None``, look for
        all l values with npowers or more modes associated with
        them

      :param l_targets: input list of l values
      :type l_targets: list of int

      :param npowers: the number of powers in the fit
      :type npowers: int

      :return: new list of l values with the above properties
      :rtype: list of int
      """

      # copy l_targets list to avoid modifying argument.  Make sure
      # each l value only occurs once in the copy.  Also, handle
      # the case l_targets is None:
      l_list = []
      if (l_targets is None):
          for mode in self.modes:
              if (mode.l not in l_list): l_list.append(mode.l)
      else:
          for l in l_targets:
              if (l not in l_list): l_list.append(l)
      
      # count the number of modes with different l values
      l_num = [0.0,]*len(l_list)
      for mode in self.modes:
          if (mode.l in l_list): l_num[l_list.index(mode.l)] += 1
      
      # remove l values from list if there are less than 2 modes with
      # that value (iterate backwards to avoid problems with shifting
      # indices).
      for i in range(len(l_list)-1,-1,-1):
          if (l_num[i] < npowers): del l_num[i], l_list[i]

      l_list.sort()
      return l_list

  def dot_product_whosglad(self, v1, v2):
      """
      Dot product used for constructing an orthonormal set of seismic indicators
      using the WhoSGlAd method (Farnir et al. 2019).

      :param v1: a vector of frequencies (or some polynomial of the radial order)
      :type v1: float np.array

      :param v2: a vector of frequencies (or some polynomial of the radial order)
      :type v2: float np.array

      :return: the dot product between v1 and v2
      :rtype: float
      """

      result = 0.0
      for i in range(len(self.modes)):
          result += v1[i]*v2[i]/self.modes[i].sigma**2
      return result

  def construct_whosglad_basis(self,maxpower=2,glitchpower=-5):
      """
      Construct orthogonal basis for the smooth component of a pulsation spectrum
      using the Gram-Schmidt orthonormalisation method as described in the WhoSGlAd
      method (Farnir et al. 2019).

      :param maxpower: maximum power in Pj(n)=n^j polynomial
      :type maxpower: int

      :param glitchpower: lower power used in Helium glitch component
      :type glitchpower: int
      """

      # easy exit:
      if ((self.Rkkinv is not None) and (self.qk is not None)): return

      l_list = self.find_l_list(None,npowers=maxpower+1)
      n = len(self.modes)
      nl = len(l_list)
      nt = (maxpower+1)*nl + 6
      pk = np.zeros((n,nt),dtype=np.float64)
      j = 0
      for l in l_list:
          for k in range(maxpower+1):
              for i in range(n):
                  if (self.modes[i].l != l): continue
                  pk[i,j] = self.modes[i].n**k
              j += 1

      j = (maxpower+1)*nl
      for i in range(n):
          ntilde = self.modes[i].n + self.modes[i].l/2.0
          pk[i,j]   = math.sin(4.0*math.pi*ntilde*config.T_He)*ntilde**glitchpower
          pk[i,j+1] = math.cos(4.0*math.pi*ntilde*config.T_He)*ntilde**glitchpower
          pk[i,j+2] = math.sin(4.0*math.pi*ntilde*config.T_He)*ntilde**(glitchpower+1)
          pk[i,j+3] = math.cos(4.0*math.pi*ntilde*config.T_He)*ntilde**(glitchpower+1)
          pk[i,j+4] = math.sin(4.0*math.pi*ntilde*config.T_CZ)*ntilde**-2.0
          pk[i,j+5] = math.cos(4.0*math.pi*ntilde*config.T_CZ)*ntilde**-2.0

      Rkk = np.zeros((nt,nt),dtype=np.float64)
      self.Rkkinv = np.zeros((nt,nt),dtype=np.float64)
      self.qk = np.zeros((n,nt),dtype=np.float64)
      for i in range(nt):
          self.qk[:,i] = pk[:,i].copy()
          for j in range(i):
              Rkk[j,i] = self.dot_product_whosglad(pk[:,i],self.qk[:,j])
              self.qk[:,i] -= self.qk[:,j]*Rkk[j,i]
          aux = math.sqrt(self.dot_product_whosglad(self.qk[:,i],self.qk[:,i]))
          self.qk[:,i] /= aux
          Rkk[i,i] = self.dot_product_whosglad(self.qk[:,i],pk[:,i])
          self.Rkkinv[i,i] = 1.0/Rkk[i,i]
          for j in range(i):
              for k in range(j,i):
                  self.Rkkinv[j,i] -= Rkk[k,i]*self.Rkkinv[j,k]/Rkk[i,i]

      # sanity check:
      freq_fit = np.zeros((n,),dtype=np.float64)
      freq_vec = np.array([mode.value for mode in self.modes])
      for j in range(nt):
          akl = self.dot_product_whosglad(self.qk[:,j],freq_vec)
          freq_fit += akl*self.qk[:,j]
          if (j == (maxpower+1)*nl-1):
              dfreq1 = freq_fit-np.array(freq_vec)
          if (j == (maxpower+1)*nl+3):
              dfreq2 = freq_fit-np.array(freq_vec)
      dfreq3 = freq_fit-np.array(freq_vec)
      self.fit = [Mode(self.modes[i].n,self.modes[i].l,f,self.modes[i].sigma) for i,f in enumerate(freq_fit)]
      chis = self.dot_product_whosglad(dfreq1,dfreq1)
      chiHe = self.dot_product_whosglad(dfreq2,dfreq2)
      chiHeCZ = self.dot_product_whosglad(dfreq3,dfreq3)
      print("{:>12}: {:e}".format("chi2(smooth)",chis))
      print("{:>12}: {:e}".format("chi2(He)",chiHe))
      print("{:>12}: {:e}".format("chi2(He+CZ)",chiHeCZ))
      return chis,chiHe,chiHeCZ

  def compute_whosglad_dnu_constraint(self,l_targets=None, maxpower=2):
      """
      Computes the average large frequency separation for a set of
      spherical degrees <l_targets>.

      :param l_targets: specifies for which l values the large frequency
             separation is to be calculated.  If ``None`` is supplied, 
             modes will be used to build an averaged value.
      :type l_targets: list of int

      :param maxpower: maximum power in Pj(n)=n^j polynomial
      :type maxpower: int
      """

      # If the basis does not exist, construct it
      if ((self.Rkkinv is None) and (self.qkinv is None)): self.construct_whosglad_basis()

      l_list_ref = self.find_l_list(None,npowers=maxpower+1)
      l_list = self.find_l_list(l_targets,npowers=maxpower+1)

      # easy exit:
      if (len(l_list) == 0):
          print("WARNING: unable to compute WhoSGlAd large separation constraint.")
          print("         Try computeing observed frequencies.")
          return

      n = len(self.modes)
      dnu = 0
      dnuSig = 0

      if (len(l_list) == 1):
        name = 'Dnu{:d}'.format(l_list[0])
        ndx = (maxpower+1)*l_list_ref.index(l_list[0])+1
        for i in range(n):
          dnu += self.qk[i,ndx]*self.Rkkinv[ndx,ndx]/(self.modes[i].sigma**2)*self.modes[i].value
        dnuSig = self.Rkkinv[ndx,ndx]
      else:
        indices = [(maxpower+1)*l_list_ref.index(l)+1 for l in l_list]
        name = 'Dnu'

        den = 0.0 
        for ndx in indices:
          den += 1.0/self.Rkkinv[ndx,ndx]**2
        for i in range(n):
          num = 0.0
          for ndx in indices:
            num += self.qk[i,ndx]/(self.Rkkinv[ndx,ndx]*self.modes[i].sigma**2)
          dnu += num/den*self.modes[i].value
        dnuSig = 1.0/np.sqrt(den)

      self.indicators[name] = Indicator(dnu,dnuSig)

  def compute_whosglad_ratio_constraint(self, l_target, maxpower=2):
      """
      Computes the normalised small separation (mean of local value)
      at a given spherical degree <l_target>.
         
      :param l_target: specifies for which l to calculate the average
             frequency ratio.
      :type l_target: int

      :param maxpower: maximum power in Pj(n)=n^j polynomial
      :type maxpower: int
      """

      l_list = self.find_l_list(None,npowers=maxpower+1)

      # easy exit:
      if (l_target == 0):
          print("WARNING: unable to compute WhoSGlAd ratio constraint using l=0.")
          print("         Please try a different constraint ...")
          return
      if ((l_target not in l_list) or (0 not in l_list)):
          print("WARNING: unable to compute WhoSGlAd ratio constraint.")
          print("         Try computeing observed frequencies.")
          return

      n = len(self.modes)
      ndx0 = (maxpower+1)*l_list.index(0)
      ndxl = (maxpower+1)*l_list.index(l_target)
      
      freq_vec = np.array([mode.value for mode in self.modes])
      a00 = self.dot_product_whosglad(self.qk[:,ndx0],freq_vec)
      a01 = self.dot_product_whosglad(self.qk[:,ndx0+1],freq_vec)
      al0 = self.dot_product_whosglad(self.qk[:,ndxl],freq_vec)
      sig = (self.Rkkinv[ndx0,ndx0]**2+self.Rkkinv[ndxl,ndxl]**2) \
          / ((a01*self.Rkkinv[ndx0+1,ndx0+1])**2) \
          + self.Rkkinv[ndx0+1,ndx0+1]**2*((a00*self.Rkkinv[ndx0,ndx0] \
          - al0*self.Rkkinv[ndxl,ndxl]) \
          / (a01*self.Rkkinv[ndx0+1,ndx0+1])**2)**2

      num = sum([self.modes[i].value*(self.qk[i,ndx0]*self.Rkkinv[ndx0,ndx0] \
              -self.qk[i,ndxl]*self.Rkkinv[ndxl,ndxl])/ \
              (self.modes[i].sigma**2) for i in range(n)])
      den= sum([self.modes[i].value*self.qk[i,ndx0+1] \
             * self.Rkkinv[ndx0+1,ndx0+1]/(self.modes[i].sigma**2) \
              for i in range(n)])
      offset = l_target/2.0 - self.Rkkinv[ndxl,ndxl+1]/self.Rkkinv[ndxl+1,ndxl+1] + self.Rkkinv[ndx0,ndx0+1]/self.Rkkinv[ndx0+1,ndx0+1]
      self.indicators['r0{:d}'.format(l_target)] = Indicator(num/den+offset,np.sqrt(sig))

  def compute_whosglad_Delta_constraint(self, l_target, maxpower=2):
      """
      Computes the ratio of large frequency separations (slope of local
      small separation) at a given spherical degree <l_target>
         
      :param l_target: specifies for which l to calculate the average
             Delta constraint.
      :type l_target: int

      :param maxpower: maximum power in Pj(n)=n^j polynomial
      :type maxpower: int
      """

      l_list = self.find_l_list(None,npowers=maxpower+1)

      # easy exit:
      if (l_target == 0):
          print("WARNING: unable to compute WhoSGlAd Delta constraint using l=0.")
          print("         Please try a different constraint ...")
          return
      if ((l_target not in l_list) or (0 not in l_list)):
          print("WARNING: unable to compute WhoSGlAd Delta constraint.")
          print("         Try computeing observed frequencies.")
          return

      sig = 0
      dnul = 'Dnu{:d}'.format(l_target)
      if self.indicators['Dnu0'].value == None: self.compute_whosglad_dnu_constraint()
      if self.indicators[dnul].value == None: self.compute_whosglad_dnu_constraint(l_targets=[l_target])
      sig = self.indicators['Dnu0'].sigma**2 \
          *(self.indicators[dnul].value/self.indicators['Dnu0'].value**2)**2\
          + self.indicators[dnul].sigma**2/self.indicators['Dnu0'].value**2
      self.indicators['Delta0{:d}'.format(l_target)] = Indicator(self.indicators[dnul].value/self.indicators['Dnu0'].value-1,np.sqrt(sig))

  def compute_whosglad_eps_constraint(self, l_target=None, maxpower=2):
      """
      Computes the epsilon (independent term) at a given spherical
      degree <l_target>.
         
      :param l_target: specifies for which l to calculate the average
             eps constraint. When None, compute the mean over all l.
      :type l_target: int

      :param maxpower: maximum power in Pj(n)=n^j polynomial
      :type maxpower: int
      """

      l_list = self.find_l_list(None,npowers=maxpower+1)
      freq_vec = np.array([mode.value for mode in self.modes])
      if l_target != None:
        # easy exit:
        if ((l_target not in l_list) or (0 not in l_list)):
            print("WARNING: unable to compute WhoSGlAd eps constraint.")
            print("         Try computeing observed frequencies.")
            return

        key = 'eps{:d}'.format(l_target)
        ndxl = (maxpower+1)*l_list.index(l_target)
        
        al0 = self.dot_product_whosglad(self.qk[:,ndxl],freq_vec)
        al1 = self.dot_product_whosglad(self.qk[:,ndxl+1],freq_vec)
        value = (al0*self.Rkkinv[ndxl,ndxl]+al1*self.Rkkinv[ndxl,ndxl+1])\
              / al1/self.Rkkinv[ndxl+1,ndxl+1]-0.5*l_target
        sig = (self.Rkkinv[ndxl,ndxl]**2+self.Rkkinv[ndxl,ndxl+1]**2\
               +(al0*self.Rkkinv[ndxl,ndxl]+al1*self.Rkkinv[ndxl,ndxl+1])**2\
                /al1**2)/(al1*self.Rkkinv[ndxl+1,ndxl+1])**2
      else:
        key = 'eps'
        value = 0
        num = 0 # num/den = mean nu
        den = 0
        sig = 0
        ndxl = lambda l: (maxpower+1)*l_list.index(l)
        ali = lambda l,i: self.dot_product_whosglad(self.qk[:,ndxl(l)+i],freq_vec)
        # num/den = mean nu 
        num = sum([ali(l,0)/self.Rkkinv[ndxl(l),ndxl(l)] for l in l_list])
        den = sum([1./self.Rkkinv[ndxl(l),ndxl(l)]**2 for l in l_list])
        numean = num/den
        # mean n+l/2
        nt = sum([0.5*l/self.Rkkinv[ndxl(l),ndxl(l)]**2 \
                  -self.Rkkinv[ndxl(l),ndxl(l)+1] \
                  /(self.Rkkinv[ndxl(l)+1,ndxl(l)+1] \
                  *self.Rkkinv[ndxl(l),ndxl(l)]**2) for l in l_list])
        nt /= den
        chin = sum([(l/(2.*self.Rkkinv[ndxl(l),ndxl(l)]))**2 \
                    -l*self.Rkkinv[ndxl(l),ndxl(l)+1] \
                    /(self.Rkkinv[ndxl(l)+1,ndxl(l)+1] \
                    *self.Rkkinv[ndxl(l),ndxl(l)]**2) \
                    +(self.Rkkinv[ndxl(l)+1,ndxl(l)+2] \
                    *self.Rkkinv[ndxl(l),ndxl(l)+1] \
                    /self.Rkkinv[ndxl(l)+1,ndxl(l)+1] \
                    -self.Rkkinv[ndxl(l),ndxl(l)+2]) \
                    /(self.Rkkinv[ndxl(l)+2,ndxl(l)+2] \
                    *self.Rkkinv[ndxl(l),ndxl(l)]**2) \
                    - nt**2/self.Rkkinv[ndxl(l),ndxl(l)]**2 for l in l_list])
        nunt = sum([ali(l,1)/self.Rkkinv[ndxl(l)+1,ndxl(l)+1] \
                    +0.5*l*ali(l,0)/self.Rkkinv[ndxl(l),ndxl(l)] \
                    -ali(l,0)*self.Rkkinv[ndxl(l),ndxl(l)+1] \
                    /self.Rkkinv[ndxl(l),ndxl(l)] \
                    /self.Rkkinv[ndxl(l)+1,ndxl(l)+1] for l in l_list])
        # Averaged dnu
        Dmean = (nunt-nt*numean*sum([1./self.Rkkinv[ndxl(l),ndxl(l)]**2 \
                 for l in l_list]))/chin
        value = numean/Dmean-nt
        sDmean = sum([((l/2.)**2-l*self.Rkkinv[ndxl(l),ndxl(l)+1] \
                      /self.Rkkinv[ndxl(l)+1,ndxl(l)+1] \
                      +(self.Rkkinv[ndxl(l)+1,ndxl(l)+2] \
                        *self.Rkkinv[ndxl(l),ndxl(l)+1] \
                        /self.Rkkinv[ndxl(l)+1,ndxl(l)+1] \
                       -self.Rkkinv[ndxl(l),ndxl(l)+2]) \
                        / self.Rkkinv[ndxl(l)+2,ndxl(l)+2] \
                       - nt**2)/(self.Rkkinv[ndxl(l),ndxl(l)]**2) \
                       for l in l_list])/(chin**2)
        sig = 1./(den*Dmean**2)+sDmean*numean**2/(Dmean**4)

      self.indicators[key] = Indicator(value,np.sqrt(sig))

  def compute_whosglad_AHe_constraint(self, maxpower=2):
      """
      Computes the helium glitch amplitude.
         
      :param maxpower: maximum power in Pj(n)=n^j polynomial
      :type maxpower: int
      """

      l_list = self.find_l_list(None,npowers=maxpower+1)
      nl = len(l_list)
      ndx = (maxpower+1)*nl # Starting index for He basis elements
      freq_vec = np.array([mode.value for mode in self.modes])
      val = sum([(self.dot_product_whosglad(self.qk[:,ndx+i],freq_vec))\
                 **2 for i in range(4)])
      self.indicators['AHe'] = Indicator(np.sqrt(val),1.)

  def compute_whosglad_ACZ_constraint(self, maxpower=2):
      """
      Computes the base of the convection zone glitch amplitude.
         
      :param maxpower: maximum power in Pj(n)=n^j polynomial
      :type maxpower: int
      """

      l_list = self.find_l_list(None,npowers=maxpower+1)
      nl = len(l_list)
      ndx = (maxpower+1)*nl+4 # Starting index for He basis elements
      freq_vec = np.array([mode.value for mode in self.modes])
      val = sum([(self.dot_product_whosglad(self.qk[:,ndx+i],freq_vec))\
                 **2 for i in range(2)])
      self.indicators['ACZ'] = Indicator(np.sqrt(val),1.)

  def compute_indicators(self):
    """
    Computes the complete set of defined seismic indicators. For an 
    understanding of their interpretation, see Farnir et al. 2019
    """
    l_list = self.find_l_list(None)
    self.compute_whosglad_dnu_constraint(l_targets=l_list)
    nl = len(l_list)
    for l in l_list:
      self.compute_whosglad_dnu_constraint(l_targets=[int(l)])
      if l!=0: 
        self.compute_whosglad_ratio_constraint(l_target=int(l))
        self.compute_whosglad_Delta_constraint(l_target=int(l))
      self.compute_whosglad_eps_constraint(l_target=int(l))
    self.compute_whosglad_eps_constraint()
    self.compute_whosglad_AHe_constraint()
    self.compute_whosglad_ACZ_constraint()

  def print_indicators(self,prefix='results',save=True):
    """
    Prints the computed indicators

    :param prefix: base name to save results as
    :type prefix: str

    :param save: Whether to save the indicators or not. It will be saved
                 as <prefix>+'-indicators.txt'
    :type save: bool
    """
    print('Computed indicators:')
    lines = ['#Indicators computed from '+prefix+'\n']
    for key in self.indicators:
      if self.indicators[key].value == None: continue # Skip undefined indicators
      line = '{:>8s}: {:7.4e} +- {:7.4e}'.format(key,self.indicators[key].value,self.indicators[key].sigma)
      print(line)
      lines = np.append(lines,line+'\n')
    if save:
      with open(prefix+'-indicators.txt','w') as f:
        f.writelines(lines)

  def echelle(self,l_targets=None,colors=pltu.colors,prefix='results',save=True):
    """"
    Plot the echelle diagram
       
    :param l_targets: list of l values to be plotted
    :type l: int

    :param colors: list of colors that will be used for plotting
    :type colors: str

    :param prefix: base name to save results as
    :type prefix: str

    :param save: Whether to save the plot or not. It will be saved as
                 <prefix>+'-echelle.pdf'
    :type save: bool
    """
    # Some plot parameters (not the best practice...)
    mks=8.#marker size
    lw=2.# line width

    fig,ax = plt.subplots(**pltu.fig_pars(1,1))
    l_list = self.find_l_list(l_targets)
    Dnu = self.indicators['Dnu'].value    
    if Dnu == None: self.compute_dnu_constraint()
    ncolors = len(colors)
    for i,l in enumerate(l_list):
      i %= ncolors # loop over the list of colors
      mask = [int(i) for i in range(len(self.modes)) if self.modes[i].l == l]
      nul = [self.modes[i].value for i in mask]
      fitl = [self.fit[i].value for i in mask]
      plt.scatter(nul%Dnu,nul,edgecolor=colors[i],facecolor='none',marker='o')
      plt.scatter(fitl%Dnu,fitl,edgecolor=colors[i],facecolor='none',marker='d')
      ax.set_xlabel(r'$\nu \% \Delta\nu~(\mu Hz)$')
      ax.set_ylabel(r'$\nu~(\mu Hz)$')

      legend_elements = [Line2D([0], [0], marker='o',color='w',label='Ref',markeredgecolor='k',markerfacecolor='none',markersize=mks,lw=lw)] \
      + [Line2D([0], [0], marker='d',color='w',label='Fit',markeredgecolor='k',markerfacecolor='none',markersize=mks,lw=lw)] \
      + [Line2D([0], [0], marker='o', color='w', label='l={:d}'.format(int(l)),markeredgecolor=colors[i],markerfacecolor='none', markersize=mks,lw=lw) for i,l in enumerate(l_list)]
    ax.legend(handles=legend_elements,ncol=1,loc=0,fontsize=10.)
    fig.tight_layout()
    if save:
      fig.savefig(prefix+'-echelle.'+pltu.ext,dpi=pltu.dpi)

  def glitch(self,l_targets=None,maxpower=2,glitchpower=-5,npt=100,colors=pltu.colors,prefix='results',save=True):
    """"
    Plot the isolated glitch
       
    :param l_targets: list of l values to be plotted
    :type l: int

    :param maxpower: maximum power in Pj(n)=n^j polynomial
    :type maxpower: int

    :param glitchpower: lower power used in Helium glitch component
    :type glitchpower: int

    :param npt: number of points at which the 'continuous' glitch will 
                be represented
    :param npt: int

    :param colors: list of colors that will be used for plotting
    :type colors: str

    :param prefix: base name to save results as
    :type prefix: str

    :param save: Whether to save the plot or not. It will be saved as
                 <prefix>+'-glitch.pdf'
    :type save: bool
    """

    def cont_basis(x,k,maxpower=2):
      """"
      Former basis elements as a function of a continuous variable
         
      :param x: point of evaluation
      :type x: float

      :param k: basis element index. Selects function to apply
      :type k: int

      :param maxpower: maximum power in Pj(n)=n^j polynomial
      :type maxpower: int
      """
      if k<=maxpower:
        p = x**k
      elif k>maxpower & k<=maxpower+4:
        p = (((k-maxpower)%2)*math.sin(4.0*math.pi*x*config.T_He) \
             +(1-(k-maxpower)%2)*math.cos(4.0*math.pi*x*config.T_He)) \
            *x**(-5+int((k-maxpower)/3))
      else:
        p = ((k%2)*math.sin(4.0*math.pi*x*config.T_CZ) \
             +(1-k%2)*math.cos(4.0*math.pi*x*config.T_CZ)) \
            *x**(-2)
      return p
    # Some plot parameters
    mks=8.#marker size
    lw=2.# line width

    n = len(self.modes)
    l_list = self.find_l_list(l_targets)
    nl = len(l_list)
    j = (maxpower+1)*nl

    freq_vec = np.array([mode.value for mode in self.modes])
    ndx = (maxpower+1)*nl # Starting index for He basis elements
    ai = lambda i: self.dot_product_whosglad(self.qk[:,ndx+i],freq_vec)
    ntilde = lambda i: self.modes[i].n + self.modes[i].l/2.0
    ndxl = lambda l: (maxpower+1)*l_list.index(l)
    ali = lambda l,i: self.dot_product_whosglad(self.qk[:,ndxl(l)+i],freq_vec)

    CZ = np.zeros(n) # Convection zone glitch only
    He = np.zeros(n) # Helium glitch only
    Res = np.zeros(n) # Ref. freq-smooth fit = isolated ref. glitch
    # Functions representing the 'continuous' spectrum divided in its
    # three contributions
    CZx = lambda ntilde: math.sin(4.0*math.pi*ntilde*config.T_CZ) \
       *ntilde**-2.0*(ai(4)*self.Rkkinv[ndx+4,ndx+4] \
       +ai(5)*self.Rkkinv[ndx+4,ndx+5]) \
       +math.cos(4.0*math.pi*ntilde*config.T_CZ)*ntilde**-2.0 \
       *ai(5)*self.Rkkinv[ndx+5,ndx+5]
    Hex = lambda ntilde: math.sin(4.0*math.pi*ntilde*config.T_He) \
       *sum([ai(j)*self.Rkkinv[ndx,ndx+j] for j in range(4)]) \
        *ntilde**glitchpower +math.cos(4.0*math.pi*ntilde*config.T_He) \
       *sum([ai(j)*self.Rkkinv[ndx+1,ndx+j] for j in range(1,4)]) \
        *ntilde**glitchpower +math.sin(4.0*math.pi*ntilde*config.T_He) \
       *sum([ai(j)*self.Rkkinv[ndx+2,ndx+j] for j in range(2,4)]) \
        *ntilde**(glitchpower+1)+math.cos(4.0*math.pi*ntilde*config.T_He) \
       *ai(3)*self.Rkkinv[ndx+3,ndx+3]*ntilde**(glitchpower+1)
    smoxj = lambda x,j : sum([self.Rkkinv[ndxl(0)+k,ndxl(0)+j] \
                        *cont_basis(x,k) if k<=maxpower \
                        else self.Rkkinv[ndx+k-maxpower,ndx+j-maxpower]\
                        *cont_basis(x,k) 
                        for k in range(j+1)])
    smox = lambda x : sum([smoxj(x,j)*ali(0,j) if j<=maxpower \
                           else smoxj(x,j)*ai(j-maxpower) \
                           for j in range(maxpower+6)])
    smonl = lambda n,l : sum([ali(l,i)*self.qk[n,ndxl(l)+i] for i in range(3)])
        
#    for l in l_list: # Compute the individual glitches at ref points
#      mask = [int(i) for i in range(n) if self.modes[i].l == l]
#      CZ[mask] = [CZx(ntilde(i)) for i in mask]
#      He[mask] = [Hex(ntilde(i)) for i in mask]

    nmin = min([self.modes[i].n for i in range(n)])
    nmax = max([self.modes[i].n for i in range(n)])
    step = float(nmax-nmin)/npt
    smo = np.array([smox(nmin+step*i) for i in range(npt)])
    Hes = np.array([Hex(nmin+step*i) for i in range(npt)])
    CZs = np.array([CZx(nmin+step*i) for i in range(npt)])
      
    fig,ax = plt.subplots(**pltu.fig_pars(1,1))
    ax.set_xlabel(r'$\nu~(\mu Hz)$')
    ax.set_ylabel(r'$\delta\nu~(\mu Hz)$')
    ncolors = len(colors)
    for i,l in enumerate(l_list): # Loop to plot glitch for each l value
      i %= ncolors # loop over the list of colors
      mask = [int(i) for i in range(n) if self.modes[i].l == l]
      nul = [self.modes[i].value for i in mask]
      sigl = [self.modes[i].sigma for i in mask]
      Res[mask] = [self.modes[i].value-smonl(i,l) for i in mask]
      ax.errorbar(nul,Res[mask],yerr=sigl,linestyle='none',label='l={:d}'.format(int(l)),**pltu.err_pars(c=colors[i]))
    # Continuous fitted glitch plotting
    ax.plot(smo+Hes+CZs,Hes,color='k',linestyle='--',label='He')
    ax.plot(smo+Hes+CZs,Hes+CZs,color='k',linestyle=':',label='He+CZ')
    ax.legend()
    fig.tight_layout()
    if save:
      fig.savefig(prefix+'-glitch.'+pltu.ext,dpi=pltu.dpi)

  def plot(self,l_targets=None,colors=pltu.colors,prefix='results',save=True,show=True):
    """
    Call all the plotting functions
       
    :param l_targets: list of l values to be plotted
    :type l: int

    :param colors: list of colors that will be used for plotting
    :type colors: str

    :param prefix: base name to save results as
    :type prefix: str

    :param save: Whether to save the plots or not.
    :type save: bool

    :param show: Whether to show plots
    :type show: bool
    """
    self.echelle(l_targets,colors=colors,prefix=prefix,save=save)
    self.glitch(l_targets,colors=colors,prefix=prefix,save=save)
    if show: 
      plt.show()
    plt.close()

  def save_freq(self,chi=[],l_targets=None,prefix='results'):
    """
    Saves the fitted frequencies along with the reference

    :param chi: Array containing the 3 chi^2 values of the adjustment
           (smooth, smooth+He, smooth+He+CZ)
    :type chi: float
       
    :param l_targets: list of l values to be plotted
    :type l: int
    
    :param prefix: base name to save results as
    :type prefix: str
    """
    n = len(self.modes)
    l_list = self.find_l_list(l_targets)
    lines = ['#Frequencies fitted from '+prefix+'\n'] 
    if len(chi) < 3:
      print('Not enough chi2 values provided, could not print them to file')
    else:
      chis    = chi[0]
      chiHe   = chi[1]
      chiHeCZ = chi[2]
      lines = np.append(lines,["#{:>12}: {:e}\n".format("chi2(smooth)",chis), \
             "#{:>12}: {:e}\n".format("chi2(He)",chiHe), \
             "#{:>12}: {:e}\n".format("chi2(He+CZ)",chiHeCZ)])
    lines = np.append(lines,('#'+2*'{:>4s}'+3*'{:>15s}'+'\n').format('l','n','nu_ref (muHz)','nu_fit (muHz)','sigma (mu Hz)'))
    for l in l_list: 
      mask = [int(i) for i in range(n) if self.modes[i].l == l]
      for i in mask:
        line = (2*'{:4d}'+3*'{:15.7e}'+'\n').format(int(self.modes[i].l),int(self.modes[i].n),self.modes[i].value,self.fit[i].value,self.modes[i].sigma) 
        lines = np.append(lines,line)
    with open(prefix+'-fit.txt','w') as f:
      f.writelines(lines)

def show_logo():
  """
  Prints the WhoSGlAd logo
  """
  print('                                                                ;\n'
  '                                                    _           ;;\n'
  '   __________________________________  _  ______ __| |_________ ;\';. _\n'
  ' ||             _                     | |       / _` |          ;  ;;\n'
  ' ||____________| \'__ ________________ | |   _  | (_| |_________ ; _ ;;\n'
  ' ||            |  _ \                 | |  / \  \__,_|          ;    ;;\n'
  ' ||____________| | | | ___ __________ |_| / _ \ _______________ ; __ ;;\n'
  ' || __        _|_| |_|/ _ \____   ____   / ___ \                ;   ;\'\n'
  ' ||_\ \      / /_____| (_)/ ___| / ___| /_/   \_\ _____________ ;  \'___\n'
  ' ||  \ \ /\ / /       \___\___ \| |  _                     ,;;;,;\n'
  ' ||___\ V  V /____________ ___) | |_| | __________________ ;;;;;; _____\n'
  '       \_/\_/             /____/ \____|                    `;;;;\'\n')

def __run__():
  """
  Runs the complete WhoSGlAd procedure
  """
  args = sys.argv
  show_logo()
  if len(args) < 2:
    sys.exit('Use: ./WhoSGlAd.py <freqfile>')
  freqfile = args[1]
  prefix = freqfile.split('.')[:-1][0] # remove extension
  modes = get_freq(freqfile)
  print('Modes used for fitting:')
  for mode in modes:
    mode.print_me()
  print()
  fitModes = Seismic()
  fitModes.modes = modes
  chis,chiHe,chiHeCZ = fitModes.construct_whosglad_basis()
  print()
  print('Fitted modes:')
  for fit in fitModes.fit:
    fit.print_me()
  fitModes.save_freq([chis,chiHe,chiHeCZ],l_targets=None,prefix=prefix)
  fitModes.compute_indicators()
  print()
  fitModes.print_indicators(prefix=prefix,save=True)
  fitModes.plot(l_targets=None,colors=pltu.colors,prefix=prefix,save=True)

if __name__ == '__main__':
# When used as a main script
  __run__()
