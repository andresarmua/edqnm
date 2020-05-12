#!/usr/bin/env python
import  sys
import  math
import  matplotlib.pyplot    as plt 
import  numpy as np
from    matplotlib  import animation
from    EDQNMmodule import get_parameters
from    EDQNMmodule import get_initial_time
from    EDQNMmodule import get_final_time

cascade_factor  = 0.9
filename        = sys.argv[1]

if 'd_' in filename: 
    u_file  = filename.replace('d_','u_')

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

if(t_f > final_time):
    print('final time exceeds simulation time')
    exit(1)
if(t_i < init_time):
    print('initial time smaller than initial simulation time')
    exit(1)

total_time      = final_time - init_time  
n_steps         = int(total_time/dt)
init_step       = int((t_i-init_time)/dt)
final_step      = int((t_f-init_time)/dt)
frame_interval  = int(np.ceil((final_step -init_step)/n_frames))
frames          = np.arange(init_step,final_step+1,frame_interval) 

if (check != 0):
    realk   = np.array([math.pow(2,i/F)  for i in range(k_max)])
else:
    realk   = np.array([i+1 for i in range(k_max)]) 


f = open(filename,'r')  
E_k_columns = np.loadtxt(f,usecols=1)/2
spectrums   = np.empty((0,k_max))
for i in range(n_steps):
    idx_start   = i * k_max
    idx_end     = (i+1) * k_max
    spectrum    = E_k_columns[idx_start:idx_end]
    spectrums   =  np.vstack((spectrums,spectrum))

if 'd_' in filename:
    u_f = open(u_file,'r')  
    E_k_column = np.loadtxt(u_f,usecols=1)
    final_spec = E_k_column[(n_steps-1)*k_max:n_steps*k_max]

#time = np.empty((0,n_steps))
__time = np.array([j*dt + init_time for j in range(n_steps)])
#__none = np.array([None for j in range(n_steps)])
#for i in range(n_steps):
#    _time   = np.concatenate((__time[:i],__none[i:])) 
#    time    = np.vstack((time,_time))


k_d = np.array([])
j_d = np.array([])
for i in range(n_steps):
    flag = False
    for j in range(k_max):
        Ed = spectrums[i]
        if Ed[j] > cascade_factor*final_spec[j]:
            k_d = np.append(k_d,realk[j])
            j_d = np.append(j_d,j)
            flag = True
            break        
    if flag == False:
        k_d = np.append(k_d,None)
        j_d = np.append(j_d,0)
j_d[:100] = 0 
j_d = j_d.astype(int)
mink    = realk.min()
maxk    = realk.max()
minE    = spectrums[-1].min()
maxE    = final_spec.max()
fig, [ax0,ax1]  = plt.subplots(nrows=2)
ax0.set_xlim(0.9*mink,1.1*maxk)
ax0.set_ylim(0.5*minE,2*maxE)
ax1.set_xlim(t_i,t_f)
ax1.set_ylim(0.9*mink,1.1*maxk)

ax0.grid(which = 'both',axis = 'both', color = '0.01', linestyle = '--', linewidth = 0.1)
ax1.grid(which = 'both',axis = 'both', color = '0.01', linestyle = '--', linewidth = 0.1)
ax0.set_yscale('log')
ax1.set_yscale('log',basey=2)
ax0.set_xscale('log',basex=2)
ax1.set_xscale('log')
ax0.set_xlabel('$k$')
ax1.set_xlabel('$t$')
ax0.set_ylabel('$E_d(k)$') 
ax1.set_ylabel('$k_E$') 
ax0.plot(realk,final_spec,'r',label='$E(k)$')
ax0.plot(realk,cascade_factor*final_spec,color='0.6',linestyle = ':')
x = np.linspace(t_i + 0.1*(t_f-t_i),t_f - 0.1*(t_f-t_i),100)
y = 64*np.power((x-t_i),-3/2)
ax1.plot(x, y, color='r', linestyle = ':',label='$-3/2$ scaling')
line0,  = ax0.plot([],[],lw=2)
dot,    = ax0.plot([],[],marker='o',ls='')
line1,  = ax1.plot([],[],lw=2)

line = [line0,dot,line1]
#if d_ plot spectrum of u_ during steady
#plt.plot(ar[realk, ar[1,spectrum,'r')

def init():
    line[0].set_data([],[]) 
    line[1].set_data([],[])
    line[2].set_data([],[]) 
    return line,

def animate(i):  
    t = i*dt + init_time
    fig.suptitle('Spectrum t = {:.2f}'.format(t))
    x = realk
    y = spectrums[i] 
    line[0].set_data([x,y])
    line[0].set_color('blue')
    ax0.lines[-2].set_label('$0.5 \, E_d(k)$')
    ax0.legend(loc= 'lower left')
    idx = j_d[i]
    if idx !=0:
        line[1].set_data(realk[idx],cascade_factor*final_spec[idx])
        line[1].set_color('green')
    x = __time[:i]
    y = k_d[:len(x)] 
    line[2].set_data(x,y)
    line[2].set_color('blue')
    ax1.lines[-1].set_label('$k_E(t)$')
    ax1.legend(loc= 'lower left')
    return line,      

anim=animation.FuncAnimation(fig,animate,init_func=init,frames=frames, interval=10)


# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html

#anim.save('spectrum_evolution_logtime.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()

