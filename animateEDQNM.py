#!/usr/bin/env python
import sys
import numpy as np
import math
from scipy import stats
import matplotlib.pyplot as plt 
from matplotlib import animation
from EDQNMmodule import get_parameters
from EDQNMmodule import get_initial_time
from EDQNMmodule import get_final_time
import inspect
import matplotlib

filename    = sys.argv[1]
if '_d_' in filename: 
    u_file  = filename.replace('_d_','_u_')

n_frames = 200
k_max, visc, check, F, d, dt = get_parameters(filename)
#compute number of steps

init_time   = get_initial_time(filename)
final_time  = get_final_time(filename)
try:
    t_i = float(sys.argv[2])
    t_f = float(sys.argv[3])
except:
    t_i = init_time
    t_f = final_time
#t_i         = float(raw_input('Set initial time: '))
#t_f         = float(raw_input('Set final time: '))
#frame_dt    = float(raw_input('Set time between frames(ms): '))

if(t_f > final_time):
    print('final time exceeds simulation time')
    exit(1)
if(t_i < init_time):
    print('initial time smaller than initial simulation time')
    exit(1)

total_time  = final_time - init_time  
n_steps     = int(total_time/dt)
init_step   = int((t_i-init_time)/dt)
final_step  = int((t_f-init_time)/dt)
frame_interval = int(np.ceil((final_step -init_step)/n_frames))
frames = np.arange(init_step,final_step+1,frame_interval) 


if (check != 0):
    realk   = np.array([math.pow(2,i/F)  for i in range(k_max)])
 
else:
    realk   = np.array([i+1 for i in range(k_max)]) 

f = open(filename,'r')  
E_k_columns = np.loadtxt(f,usecols=1)
spectrums   = np.empty((0,k_max))
for i in range(n_steps):
    idx_start   = i * k_max
    idx_end     = (i+1) * k_max
    spectrum    = E_k_columns[idx_start:idx_end]
    spectrums   =  np.vstack((spectrums,spectrum))

if '_d_' in filename:
    u_f = open(u_file,'r')  
    E_k_column = np.loadtxt(u_f,usecols=1)
    final_spec   = E_k_column[(n_steps-1)*k_max:n_steps*k_max]
    spectrums = spectrums/2 #otherwise Ed -> 2E, this way they match

mink    = realk.min()
maxk    = realk.max()
minE    = spectrums[-1].min()
maxE    = spectrums[-1].max()

fig     = plt.figure()
ax      = plt.axes(xlim=(0.9*mink,1.1*maxk),ylim=(0.5*minE,2*maxE))

plt.grid(which = 'both',axis = 'both', color = '0.01', linestyle = '--', linewidth = 0.1)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('$k$')
plt.ylabel('$E(k)$') if 'u_' in filename else plt.ylabel('$E_d(k)$')
#plt.title('Spectrum(t)'))
if '_d_' in filename: plt.plot(realk,final_spec,'r')
line, = ax.plot([],[],lw=2)

#if d_ plot spectrum of u_ during steady
#plt.plot(ar[realk, ar[1,spectrum,'r')
def init():

    line.set_data([],[]) 
    return line,

def animate(i):  
    
    t = i*dt + init_time
    plt.title('Spectrum t = {:.2f}'.format(t))

    x = realk
    y = spectrums[i] 
    line.set_data(x,y)
    line.set_color('blue')
    return line,      

anim=animation.FuncAnimation(fig,animate,init_func=init,frames=frames, interval=5)


# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html

#anim.save('spectrum_evolution.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()

