#!/usr/bin/env python

import numpy as np
import sys 
import matplotlib.pyplot as plt
import os 
#tf = 300
filename = sys.argv[1]
f = open(filename,'r')
data = np.loadtxt(f,usecols=(0,1))
t,E_k = data.transpose()
#t = t[:tf]
#E_k = E_k[:tf]


a = plt.figure()
plt.xlabel('t')
plt.ylabel('$E$')
plt.yscale('log')
plt.plot(t,E_k)
#plt.savefig('pert_run_Test.png')
plt.show()
