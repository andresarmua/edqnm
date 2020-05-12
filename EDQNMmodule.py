#!/usr/bin/env python

import numpy as np
from scipy.stats import linregress

def get_initial_time(ufilename):
    '''return initial time of spc file'''
    if '.stats' in ufilename:
        ufilename = ufilename.replace('stats','spc')
    with open(ufilename) as readfile:
        ulines = readfile.readlines()
    init_time = -1
    spectra_block_check = 'asd'
    for line in ulines:
        line_data = line.split()
        
        try:
            spectra_block_check = line_data[1]
        except:
            pass

        if spectra_block_check == "Plot":
            init_time = float(line_data[-1])
        
        if init_time != -1:
            break
    
    return init_time


def get_final_time(ufilename):
    '''return final time of spc file'''
    if '.stats' in ufilename:
        ufilename = ufilename.replace('stats','spc') 
    with open(ufilename) as readfile:
        ulines = readfile.readlines()
    final_time = -1
    ulines.reverse()
    spectra_block_check = 'asd'
    for line in ulines:
        line_data = line.split()
    
        try:
            spectra_block_check = line_data[1]
        except:
            pass

        if spectra_block_check == "Plot":
            final_time = float(line_data[-1])
        
        if final_time != -1:
            break
    
    return final_time




# --- get parameters from header of stats or spc files --- 

def get_parameters(ufilename):
    '''function returns lattice size(int), viscosity (float), \
logcheck (int),f (float), dimension (int),\ 
output time interval (float)'''

    with open(ufilename) as readfile:
        ulines = readfile.readlines()
    
    target_parameters   = ['KMAX','VISC','LOG','F','DIMENSION','INTERVAL']

    for line in ulines:
        
        try:
            parameter   = line.split()[1]
            value       = line.split()[2]
        except:
            print("spc file not formatted correctly in \
line {}".format(ulines.index(line)))
            exit(1)
        
        if   parameter == target_parameters[0]:
            lattice = int(value)
        elif parameter == target_parameters[1]:
            visc = float(value)
        elif parameter == target_parameters[2]:
            check = int(value)
        elif parameter == target_parameters[3]:
            F = float(value)
        elif parameter == target_parameters[4]:
            d = int(value)
        elif parameter == target_parameters[5]:
            dt = float(value)
        elif parameter == 'Plot':
            break

    return lattice,visc,check,F,d,dt




def get_viscosity(ufilename):
    '''return viscosity from header'''
    with open(ufilename) as readfile:
        ulines = readfile.readlines()
    
    target_parameter    = 'VISC'

    for line in ulines:
        
        try:
            parameter   = line.split()[1]
            value       = line.split()[2]
        except:
            print("spc file not formatted correctly in \
line {}".format(ulines.index(line)))
            exit(1)
        
        if parameter == target_parameter:
            visc = float(value)
        elif parameter == 'Plot':
            break

    return visc


def get_lattice(ufilename):
    '''return lattice (k_max) from header'''
    with open(ufilename) as readfile:
        ulines = readfile.readlines()
    
    target_parameter    = 'KMAX'

    for line in ulines:
        
        try:
            parameter   = line.split()[1]
            value       = line.split()[2]
        except:
            print("spc file not formatted correctly in \
line {}".format(ulines.index(line)))
            exit(1)
        
        if parameter == target_parameter:
            lattice = float(value)
        elif parameter == 'Plot':
            break

    return lattice



def get_logcheck(ufilename):
    '''return logcheck from header'''
    with open(ufilename) as readfile:
        ulines = readfile.readlines()
    
    target_parameter    = 'LOG'


    for line in ulines:
        
        try:
            parameter   = line.split()[1]
            value       = line.split()[2]
        except:
            print("spc file not formatted correctly in \
line {}".format(ulines.index(line)))
            exit(1)
        
        if parameter == target_parameter:
            logcheck = float(value)
        elif parameter == 'Plot':
            break

    return logcheck


def get_dimension(ufilename):
    '''return dimension from header'''
    with open(ufilename) as readfile:
        ulines = readfile.readlines()
    
    target_parameter    = 'DIMENSION'


    for line in ulines:
        
        try:
            parameter   = line.split()[1]
            value       = line.split()[2]
        except:
            print("spc file not formatted correctly in \
line {}".format(ulines.index(line)))
            exit(1)
        
        if parameter == target_parameter:
            dimension = float(value)
        elif parameter == 'Plot':
            break

    return dimension


def get_f(ufilename):
    '''return f (2^(i/f)) from header'''
    with open(ufilename) as readfile:
        ulines = readfile.readlines()
    
    target_parameter = 'F'


    for line in ulines:
        
        try:
            parameter   = line.split()[1]
            value       = line.split()[2]
        except:
            print("spc file not formatted correctly in \
line {}".format(ulines.index(line)))
            exit(1)
        
        if parameter == target_parameter:
            f = float(value)
        elif parameter == 'Plot':
            break

    return f


def get_dt(ufilename):
    '''return time intervala of the output \
(not dt of the simulation) from header'''
    with open(ufilename) as readfile:
        ulines = readfile.readlines()
    
    target_parameter = 'INTERVAL'


    for line in ulines:
        
        try:
            parameter   = line.split()[1]
            value       = line.split()[2]
        except:
            print("spc file not formatted correctly in \
line {}".format(ulines.index(line)))
            exit(1)
        
        if parameter == target_parameter:
            dt = float(value)
        elif parameter == 'Plot':
            break

    return dt




def col_L(ufilename):
    if 'd_' in ufilename:
        ufilename = ufilename.replace('d_','u_')
        print('Warning: d_ file was introduced instead of \
u_ to compute L, script fixes it if possible')
    if '.spc' in ufilename:
        ufilename = ufilename.replace('.spc','.stats')
    try:
        readfile = open(ufilename) 
    except:
        # in the future: add stats file creation from script
        print('Error: stats file missing, run script on spc file')
        exit(1)
    
    ulines = readfile.readlines()
        
    for line in ulines:
        line_data = line.split('\t')
        for data in line_data:
            if ' L ' in data: 
                idx = line_data.index(data)
    L = np.loadtxt(ufilename,usecols=idx)
        
        

    return L



def col_Re(ufilename):
    if 'd_' in ufilename:
        ufilename = ufilename.replace('d_','u_')
        print('Warning: d_ file was introduced instead of \
u_ to compute Re, script fixes it if possible')
    if '.spc' in ufilename:
        ufilename = ufilename.replace('.spc','.stats')
    try:
        readfile = open(ufilename) 
    except:
        # in the future: add stats file creation from script
        print('Error: stats file missing, run script on spc file')
        exit(1)
    
    ulines = readfile.readlines()
        
    for line in ulines:
        line_data = line.split('\t')
        for data in line_data:
            if 'Re_L ' in data: 
                idx = line_data.index(data)
    Re = np.loadtxt(ufilename,usecols=idx)
        
    return Re 

def col_U(ufilename):
    if 'd_' in ufilename:
        ufilename = ufilename.replace('d_','u_')
        print('Warning: d_ file was introduced instead of \
u_ to compute U, script fixes it if possible')
    if '.spc' in ufilename:
        ufilename = ufilename.replace('.spc','.stats')
    try:
        readfile = open(ufilename) 
    except:
        # in the future: add stats file creation from script
        print('Error: stats file missing, run script on spc file')
        exit(1)
    
    ulines = readfile.readlines()
        
    for line in ulines:
        line_data = line.split('\t')
        for data in line_data:
            if ' U ' in data: 
                idx = line_data.index(data)
    U = np.loadtxt(ufilename,usecols=idx)
        
    return U



def col_diss_rate(ufilename):
    if 'd_' in ufilename:
        ufilename = ufilename.replace('d_','u_')
        print('Warning: d_ file was introduced instead of \
u_ to compute diss, script fixes it if possible')
    if '.spc' in ufilename:
        ufilename = ufilename.replace('.spc','.stats')
    try:
        readfile = open(ufilename) 
    except:
        # in the future: add stats file creation from script
        print('Error: stats file missing, run script on spc file')
        exit(1)
    
    ulines = readfile.readlines()
        
    for line in ulines:
        line_data = line.split('\t')
        for data in line_data:
            if ' diss ' in data: 
                idx = line_data.index(data)
    diss = np.loadtxt(ufilename,usecols=idx)
        
    return diss


def col_t(ufilename):
    if '.spc' in ufilename:
        ufilename = ufilename.replace('.spc','.stats')
    try:
        readfile = open(ufilename) 
    except:
        # in the future: add stats file creation from script
        print('Error: stats file missing, run script on spc file')
        exit(1)
    
    ulines = readfile.readlines()
        
    for line in ulines:
        line_data = line.split('\t')
        for data in line_data:
            if ' time ' in data: 
                idx = line_data.index(data)
    time = np.loadtxt(ufilename,usecols=idx)
        
    return time

def col_E(ufilename):
    if '.spc' in ufilename:
        ufilename = ufilename.replace('.spc','.stats')
    try:
        readfile = open(ufilename) 
    except:
        # in the future: add stats file creation from script
        print('Error: stats file missing, run script on spc file')
        exit(1)
    
    ulines = readfile.readlines()
        
    for line in ulines:
        line_data = line.split('\t')
        for data in line_data:
            if ' energy ' in data: 
                idx = line_data.index(data)
    energy = np.loadtxt(ufilename,usecols=idx)
        
    return energy

# --- compute mean and std of parameters in the steady state ---

def get_L(ufilename):
    '''returns the mean and std value of L using the \
last 1000 values of the L column in u stats file'''
    #has to improve to avoid hard coding the 1000
    #one option is to use the std to guess when the 
    #steady state has started
    _L = col_L(ufilename)
    _L = _L[-1000:]
    L = _L.mean()
    std_L= _L.std()
    return L, std_L

def get_E(ufilename):
    '''returns the mean and std value of L using the \
last 1000 values of the L column in u stats file'''
    #has to improve to avoid hard coding the 1000
    #one option is to use the std to guess when the 
    #steady state has started
    _E = col_E(ufilename)
    _E = _E[-1000:]
    E = _E.mean()
    std_E= _E.std()
    return E, std_E

def get_Re(ufilename):
    '''returns the mean and std value of Re using the \
last 1000 values of the Re_L column in u stats file'''
    #has to improve to avoid hard coding the 1000
    #one option is to use the std to guess when the 
    #steady state has started
    _Re = col_Re(ufilename)
    _Re = _Re[-1000:]
    Re = _Re.mean()
    std_Re= _Re.std()
    return Re, std_Re


def get_U(ufilename):
    '''returns the mean and std value of U_rms using the \
last 1000 values of the U column in u stats file'''
    #has to improve to avoid hard coding the 1000
    #one option is to use the std to guess when the 
    #steady state has started
    _U = col_U(ufilename)
    _U = _U[-1000:]
    U = _U.mean()
    std_U = _U.std()
    return U, std_U


def get_diss_rate(ufilename):
    '''returns the mean and std value of dissipation rate \
using the last 1000 values of the diss column in u stats file'''
    #has to improve to avoid hard coding the 1000
    #one option is to use the std to guess when the 
    #steady state has started
    _diss = col_diss_rate(ufilename)
    _diss = _diss[-1000:]
    diss = _diss.mean()
    std_diss= _diss.std()
    return diss, std_diss


def get_T(ufilename):
    '''returns the mean and std value of the large eddy turnover \
            time using the values of U and L'''
    U, stdU = get_U(ufilename)
    L, stdL = get_L(ufilename)
    T       = L/U
    std_T   = T*(stdU / U + stdL / L)
    return T, std_T

def get_Te(ufilename):
    '''returns the mean and std value of the large eddy turnover \
time using the values of E and dissipation rate'''
    #steady state has started
    E, stdE = get_E(ufilename)
    eps, stdeps = get_diss_rate(ufilename)
    Te      = E/eps
    std_Te  = Te*(stdE / E + stdeps / eps)
    return Te, std_Te

def get_tau(ufilename):
    '''returns the mean and std value of the Kolmogorov  \
time microscale using the viscosity and dissipation rate'''
    nu = get_viscosity(ufilename)
    eps, stdeps = get_diss_rate(ufilename)
    tau     = (nu/eps)**0.5
    std_tau   = stdeps * nu**0.5 / (2*eps**1.5)

    return tau, std_tau

def get_eta(ufilename):
    '''returns the mean and std value of the Kolmogorov  \
length microscale using the viscosity and dissipation rate'''
    nu = get_viscosity(ufilename)
    eps, stdeps = get_diss_rate(ufilename)
    eta         = (nu**3/eps)**0.25
    std_eta     = stdeps * np.power(nu,3/4) / (4*np.power(eps,5/4))

    return eta, std_eta

def get_lyapunov(ufilename):
    '''returns the mean and std value of lyapunov exponent''' 
    if 'u_' in ufilename:
        ufilename = ufilename.replace('u_','d_')
        print('Warning: u_ file was introduced instead of d_')
    t = col_t(ufilename) 
    E = col_E(ufilename)
    t_i,t_f,i_i,i_f = get_limits(t,E)
#    t = t[i_i:i_f]
#    E = E[i_i:i_f]  
    t = t[300:700]
    E = E[300:700]
    try:
        y = np.log(E)
    except:
        print('non positive energy in {}'.format(ufilename))
        exit(1)
    X = t

    _lyapunov, intercept, rvalue, pvalue, _std = linregress(X,y)
    lyapunov    = _lyapunov/2 #because we used E, not u
    std_lyap    = _std/2  
    return lyapunov,std_lyap







def sec_d(t,E):  
    lE = np.log10(E)
    second_derivative = np.array([])
    for i in range(len(t)):
        if i == 0:
            _second = (lE[i+2] - 2 * lE[i+1] + lE[i])/(t[i+1]-t[i])**2
        elif (i!=0 and i!=len(t)-1):
            _second = (lE[i+1] - 2 * lE[i] + lE[i-1])/(t[i]-t[i-1])**2
        elif (i == len(t)-1):
            _second = (lE[i] - 2 * lE[i-1] + lE[i-2])/(t[i]-t[i-1])**2
        second_derivative = np.append(second_derivative,_second)
    return second_derivative


def fir_d(t,E):   
    lE = np.log10(E)
    first_derivative = np.array([])
    for i in range(len(t)):
        if i == 0:
            _first = (lE[i+1] - lE[i])/(t[i+1]-t[i])
        else:
            _first = (lE[i] - lE[i-1])/(t[i]-t[i-1])
        first_derivative = np.append(first_derivative,_first)
    return first_derivative



def get_limits(t,E):
    #hard_coding by looking at the plots
    #and considering dt(interval) = 0.01
    second = sec_d(t,E)
    flag_i = True
    t_i = t[100]
    i_i = 100
    i_f = len(t)
    t_f = t[len(t)-1]
    for i in range(100,len(t)):
        if (second[i] < 0.05 and flag_i):
            i_i = i
            t_i = t[i]
            flag_i = False
        elif second[i] < -0.05:
            t_f = t[i]
            i_f = i
            if (t_f - t_i) < 0.3:
                t_f = t_f + 0.5
                i_f = i_f + 50 
            break
    return t_i, t_f,i_i,i_f
