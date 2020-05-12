#!/usr/bin/python
import sys
import numpy as np
import math
from scipy import stats
import matplotlib.pyplot as plt 
from matplotlib import animation
from EDQNMmodule import get_parameters
from EDQNMmodule import get_initial_time
from EDQNMmodule import get_final_time


filename    = sys.argv[1]


k_max, visc, check, F, d, dt = get_parameters(filename)
#compute number of steps

init_time   = get_initial_time(filename)
final_time  = get_final_time(filename)

try:
    t = float(sys.argv[2])
except:
    t = init_time
#t_i         = float(raw_input('Set initial time: '))
#t_f         = float(raw_input('Set final time: '))
#frame_dt    = float(raw_input('Set time between frames(ms): '))

if(t > final_time or t < init_time):
    print('time not present in simulation')
    exit(1)

total_time  = final_time - init_time  
n_steps = int(total_time/dt)


i = int(t/dt)
if (check != 0):
    realk   = np.array([math.pow(2,i/F)  for i in range(k_max)])
 
else:
    realk   = np.array([i+1 for i in range(k_max)]) 

f = open(filename,'r')  

E_k_column  = np.loadtxt(f,usecols=1)
spectrums   = np.empty((0,k_max))
for i in range(n_steps):
    idx_start   = i * k_max
    idx_end     = (i+1) * k_max
    spectrum    = E_k_column[idx_start:idx_end]
    spectrums   =  np.vstack((spectrums,spectrum))


mink    = realk.min()
maxk    = realk.max()
minE    = spectrums[i].clip(1e-33).min()
maxE    = spectrums[i].max()


fig     = plt.figure()
ax      = plt.axes(xlim=(0.9*mink,1.1*maxk),ylim=(0.9*minE,1.1*maxE))

plt.grid(which = 'both',axis = 'both', color = '0.01', linestyle = '--', linewidth = 0.1)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('$k$')
plt.ylabel('$E(k)$') if 'u_' in filename else plt.ylabel('$E_d(k)$')
#plt.title('Spectrum(t)'))

#if d_ plot spectrum of u_ during steady
#plt.plot(ar[realk, ar[1,spectrum,'r')

    
t = i * dt
plt.title('Spectrum t = %.2f' %t)

x = realk
y = spectrums[0] 
plt.plot(x,y,'.')

plt.show()

