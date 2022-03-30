#!/usr/bin/env python
# This module defines useful plotting utilities to keep constant conventions
# such as font size and color schemes.
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D

plt.rc('text',usetex=True)
plt.rc('font',family='serif')
plt.rc('legend',fontsize=10.5)
plt.rc('xtick',labelsize=18)
plt.rc('ytick',labelsize=18)
plt.rc('axes',labelsize=20)
plt.rc('axes',titlesize=20)
plt.rc('savefig',dpi=300)
plt.rc('savefig',format='pdf')

def make_error_boxes(ax,x,y,xerr,yerr,alpha=1.,fill=False,ec='k',fc='k',lw=1.,label='errbox'):
  import mimic_alpha as malp
  import numpy as np
  xmin, xmax = ax.get_xlim()
  ymin, ymax = ax.get_ylim()
  ax.add_patch(
    patches.Rectangle(
#        ((x-xerr-xmin)/(xmax-xmin),(y-yerr-ymin)/(ymax-ymin)),
#        (2*xerr)/(xmax-xmin),
#        (2*yerr)/(ymax-ymin),
        ((x-xerr),(y-yerr)),
        (2*xerr),
        (2*yerr),
        lw=lw,
        fill=fill,      # background
        alpha=alpha,    # transparancy
        edgecolor=ec,   # borders
        facecolor=np.ravel(malp.colorAlpha_to_rgb(fc,alpha, bg='w'))    # background
    )
  )
  ax.scatter(x,y,marker='+',color=fc,lw=1.,label=label)
  return ax

def slanted_error_boxes(ax,x,R,xerr,Rerr,alpha=1.,fill=False,ec='black',fc='b',lw=1.,label='Obs'):
  import mimic_alpha as malp
  import numpy as np
  from matplotlib.path import Path
  from matplotlib.patches import PathPatch
  def func_edges(R,x):
    y = 4.*np.log10(x/5777.)+2.*np.log10(R)
#    print(y,R,np.log10(R),x,np.log10(x/5777.))
    return y
  xmin, xmax = ax.get_xlim()
  ymin, ymax = ax.get_ylim()
  codes = [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY]
  vertices = [(np.log10(x+xerr), func_edges(R-Rerr,x+xerr)), (np.log10(x-xerr), func_edges(R-Rerr,x-xerr)), (np.log10(x-xerr), func_edges(R+Rerr,x-xerr)), (np.log10(x+xerr), func_edges(R+Rerr,x+xerr)), (np.log10(x+xerr), func_edges(R-Rerr,x+xerr))]
  vertices = np.array(vertices, float)
  path = Path(vertices, codes)
  pathpatch = PathPatch(path, facecolor='None', edgecolor=ec)
  ax.add_patch(pathpatch)
  ax.scatter(np.log10(x),func_edges(R,x),marker='+',color=fc,lw=1.,label=label)
  return ax

def make_err_range(ax,val,err,c='k',w='v',label='Obs'):
  if w=='v':
    ax.axvline(val,0,1,color=c,label=label)
    ax.axvline(val+err,0,1,color=c,linestyle='--')
    ax.axvline(val-err,0,1,color=c,linestyle='--')
  elif w=='h':
    ax.axhline(val,0,1,color=c,label=label)
    ax.axhline(val+err,0,1,color=c,linestyle='--')
    ax.axhline(val-err,0,1,color=c,linestyle='--')
  else:
    print('Unknown keyword'+w+' for w parameter. Should either be v or h.')
  return ax

def fig_pars(nrows=1,ncols=1):
  return {'nrows':nrows,'ncols':ncols,'figsize':(ncols*6,nrows*4)}

def err_pars(c='k',mk=None,lw=.5):
  return {'c':c,'fmt':'none','marker':mk,'lw':lw,'zorder':0,'capsize':2}

label_kwargs={'fontsize':14}
ticks_kwargs={'axis':'both','which':'major','labelsize':12}
cbar_kwargs={'labelpad':-20,'y':1.1,'rotation':0,'fontsize':14}
ticklabel_kwargs={'format':'sci','scilimits':(-3,4),'axis':'both'}
ext='pdf'
dpi=200

colors = ('#204473','#4393D9','#BC7C2A','#F2D1B3','#592512','#A4B3BF','#F2A2B1','#6732D9','#F2C849','#F26D3D','#D95284','#025959','#BF2C1F','#59522F','#A69958','#BFA29B','#ACF2DB')
mk_list = ('d','.','o','v','^','<','>','s','p')
