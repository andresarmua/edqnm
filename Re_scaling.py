#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import EDQNMmodule as ed
import sys
import os
import pandas as pd
import statsmodels.api as sm
from glob import glob
from scipy.stats import linregress
from   cycler import cycler

def main(): 
    try:
        directory = sys.argv[1]
    except:
        directory = ''
    filenames   =   '*u_*.stats'
    Path = os.path.join(directory,filenames)
    files       = glob(Path)

    n_files     = len(files)

    Re          = np.array([])
    lyapunov    = np.array([])
    T           = np.array([])
    Te          = np.array([])
    viscosity   = np.array([])
    lattice     = np.array([])
    F           = np.array([])
    check       = np.array([])
    k_max       = np.array([])
    eps         = np.array([])
    L           = np.array([])
    tau         = np.array([])
    eta         = np.array([])

    std_Re      = np.array([])
    std_ly      = np.array([])
    std_T       = np.array([])
    std_Te      = np.array([])
    std_eps     = np.array([])
    std_L       = np.array([])
    std_tau     = np.array([])
    std_eta     = np.array([])


    for f in files:
        _Re, _std_Re    = ed.get_Re(f)
        _T,  _std_T     = ed.get_T(f)
        _Te, _std_Te    = ed.get_Te(f)
        _viscosity      = ed.get_viscosity(f)
        _lattice        = ed.get_lattice(f)
        _F              = ed.get_f(f)
        _check          = ed.get_logcheck(f) 
        _eps, _std_eps  = ed.get_diss_rate(f)
        _L,  _std_L     = ed.get_L(f)
        _tau, _std_tau  = ed.get_tau(f)
        _eta, _std_eta  = ed.get_eta(f)

        f = f.replace('u_','d_')
        _ly, _std_ly    = ed.get_lyapunov(f)
        
        
        _k_max = np.power(2,(_lattice-1)/_F) if (_check != 0) \
                else _lattice 
        
         
        Re          = np.append(Re,_Re)
        lyapunov    = np.append(lyapunov,_ly)
        T           = np.append(T,_T)
        Te          = np.append(Te,_Te)
        viscosity   = np.append(viscosity,_viscosity)
        lattice     = np.append(lattice,_lattice)
        F           = np.append(F,_F)
        check       = np.append(check,_check)
        k_max       = np.append(k_max,_k_max)
        eps         = np.append(eps,_eps)
        L           = np.append(L,_L)
        tau         = np.append(tau,_tau)
        eta         = np.append(eta,_eta)

 
        
        std_Re      = np.append(std_Re,_std_Re)
        std_ly      = np.append(std_ly,_std_ly)
        std_T       = np.append(std_T,_std_T)
        std_Te      = np.append(std_Te,_std_Te)
        std_eps     = np.append(std_eps,_std_eps)
        std_L       = np.append(std_L,_std_L)
        std_tau     = np.append(std_tau,_std_tau)
        std_eta     = np.append(std_eta,_std_eta)

    srt_idx     = Re.argsort()

    Re          = Re[srt_idx]
    lyapunov    = lyapunov[srt_idx]
    T           = T[srt_idx]
    Te          = Te[srt_idx]
    viscosity   = viscosity[srt_idx]
    lattice     = lattice[srt_idx]
    F           = F[srt_idx]
    check       = check[srt_idx]
    k_max       = k_max[srt_idx]
    eps         = eps[srt_idx]
    L           = L[srt_idx]
    tau         = tau[srt_idx]
    eta         = eta[srt_idx]

    std_Re      = std_Re[srt_idx]
    std_ly      = std_ly[srt_idx]
    std_T       = std_T[srt_idx]
    std_Te      = std_Te[srt_idx]
    std_eps     = std_eps[srt_idx]
    std_L       = std_L[srt_idx]
    std_tau     = std_tau[srt_idx]
    std_eta     = std_eta[srt_idx]

    std_lyT     = lyapunov * T * (std_ly/lyapunov + std_T/T)
    #--- create table of values
    df = pd.DataFrame( { '$N$':lattice.astype(int),'$\nu$':viscosity,\
                        'Re':Re.astype(int),'$\lambda$':lyapunov,\
                        '$\\varepsilon$':eps,'$T$':T, \
                        '$\eta k_m$': k_max*eta } )

    df = df.round({'$\nu$':5,'$\lambda$':3,'$\sigma_{\lambda}$':3,\
                    '$T$':3})
    print(df.to_latex(index = False,escape=False))
    #---

    #--- perform linear regression on scaling lyap T vs Re
    alpha,intcpt,rvalue,pvalue,std_alpha = \
            linregress(np.log10(Re), np.log10(lyapunov*T))
    print('alpha = {:5.3f} +/- {:5.3f}'.format(alpha,std_alpha)) 
    #---
    
    
    #---plot & save
    set_plot_features()
    D = 10**intcpt
    Ruelle_fit  = power_law(D,alpha)
    plot(Re, lyapunov*T, xerr = std_Re, yerr = std_lyT,\
            xscale= 'log', yscale= 'log', \
            xlabel = '$Re$', ylabel = '$\lambda T$',\
            func = Ruelle_fit,save_as = 'Ruelle_3D.png')
    #---


def set_plot_features():
    params = {  'legend.fontsize':'x-large','axes.labelsize':'xx-large',\
                'xtick.labelsize':'x-large','ytick.labelsize':'x-large'}
    plt.rcParams.update(params)
    plt.rc('text',usetex=True)
    plt.rc('font',family='serif')
#    plt.rcParams['axes.prop_cycle'] = cycler(color = 'krgcmyb')
    plt.style.use('seaborn-muted')

def linear_function(slope,intercept):
    def func(x):
       return slope * x + intercept
    return func

def power_law(A,scaling):
    def func(x):
        return A * x**scaling
    return func

def plot(x,y,xerr,yerr, func = None, func2 = None, xlabel = '', ylabel = '', \
        xscale = None, yscale = None, save_as = None ):

    fig         = plt.figure()
    if  xscale == None:
        
        x_margin    = 0.1 * (max(x) - min(x))
        x_left      = min(x) - x_margin
        x_right     = max(x) + x_margin
    
    elif xscale == 'log':
        if np.log10(max(x)/min(x)) < 1.5:
            x_left  = 0.9 * min(x)
            x_right = 1.1 * max(x)
        else:  
            x_margin    = 0.3 * (max(x) - min(x))
            x_left      = min(x) - x_margin
            x_right     = max(x) + x_margin
    if yscale == None:     
        y_margin    = 0.1 * (max(y) - min(y))
        y_bottom    = 0.*min(y) - y_margin
        y_top       = max(y) + y_margin
    elif yscale == 'log':
        if np.log10(max(y)/min(y)) < 1.5:
            y_bottom  = 0.9 * min(y)
            y_top = 1.1 * max(y)
        else:  
            y_margin    = 0.1 * (max(y) - min(y))
            y_bottom      = min(y) - y_margin
            y_top     = max(y) + y_margin
    
    
    ax          = plt.axes( xlim=(x_left, x_right),\
                    ylim=(y_bottom, y_top) )
    ax.errorbar(x, y, xerr = xerr, yerr = yerr, marker ='.', markerfacecolor = 'None',linestyle = 'None',label ='Data',ecolor = '0.6', capsize = 3, elinewidth = 1, capthick = 1)
    # ax.plot(x,y,'.') # fix the choice or not of xerr and yerr
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xscale   != None: ax.set_xscale(xscale)
    if yscale   != None: ax.set_yscale(yscale)
    x = np.append(x, [x_left, x_right])
    if func     != None: ax.plot(x,func(x),'-')
    if func2    != None: ax.plot(x,func2(x),'-')
    ax.set_xticks([400,1000,5000])
    if save_as  != None: fig.savefig(save_as, bbox_inches = 'tight', format = "png", dpi = 300)
    plt.show()

main()
