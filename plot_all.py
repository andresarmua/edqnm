#!/usr/bin/env python

import numpy as np
import sys 
import matplotlib.pyplot as plt
import os 
from glob import glob 
from EDQNMmodule import get_viscosity

try:

    u_or_d = sys.argv[1]
except:
    u_or_d = 'u'

def get_t_and_E(filename):
    f = open (filename,'r')
    data = np.loadtxt(f,usecols=(0,1))
    t,E_k = data.transpose()

    return t,E_k

try:
    directory    = sys.argv[2] 
except:
    directory = ''
filen   = '*'+u_or_d+'_*.stats'
Path    = os.path.join(directory,filen)
stats_files = glob(Path)
#tf=300


n_plots     = len(stats_files)

viscosity   = np.array([])


for f in stats_files:
    _viscosity = get_viscosity(f)
    viscosity = np.append(viscosity,_viscosity)

srt_idx         = viscosity.argsort()

fig, ax         = plt.subplots(nrows=n_plots,sharex=True)
j = 0
for i in srt_idx:
    ax[j].set_xlabel('t')
    ax[j].set_ylabel('$E$')
#    ax[j].set_ylim(bottom = 1e-8, top =1)
#    ax[j].set_xlim(left = 17.2, right = 30 )
    if u_or_d == 'd': ax[j].set_yscale('log')
    t,E_k = get_t_and_E(stats_files[i])
    ax[j].plot(t,E_k,label = 'viscosity = {:.4f}'.format(viscosity[i]))
    ax[j].legend(loc='lower right')
    j += 1
    #fig.savefig('all_E_vs_t.png')
plt.show()
