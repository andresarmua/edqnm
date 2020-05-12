#!/usr/bin/python
import math
import sys
import os
from scipy.special  import gamma


# Works like HDprocess although you need to input from the terminal the number of wavenumbers, the viscosity, if logirithmic disretisation is used and the F value.

ufilename       = sys.argv[1]
output_filename = ufilename.replace("spc","stats")

with open(ufilename) as readfile:
        ulines = readfile.readlines()

# find global variables from header, visc, viscpower, etc

target_parameters = ['KMAX','VISC','LOG','F','DIMENSION']

for line in ulines:
    
    try:
        parameter   = line.split()[1]
        value       = line.split()[2]
    except:
        print("spc file not formatted correctly in line {}".format(ulines.index(line)))
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
    elif parameter == 'Plot':
        break

if (check != 0):
    realk   = [math.pow(2,i/F)  for i in range(lattice)]
    dk      = [math.pow(2, i/F)*(math.pow(2,1.0/(2*F)) - math.pow(2,-1.0/(2*F))) for i in range(lattice)]
 
else:
    realk   = [i+1 for i in range(lattice)] 
    dk      = [1.0 for i in range(lattice)]

#create output file (delete and create if exists)
try:
    os.remove(output_filename)
except OSError:
    pass
output_file = open(output_filename,"w+")

# write the parameter list into the stats file
for line in ulines:   
    if line.split()[0] == '#':
        output_file.write(line)
    elif '#1:k' in line.split()[0]:
        break
    else:
        print('spc file not formatted correctly in line {}'.format(ulines.index(line)))
        os.remove(output_filename)
        exit(1)
kmax = lattice + 1  #actually one greater than kmax, to make loop not need - 1
factor = 2.0/float(d)

output_file.write("{:<7} \t {:<16} \t {:<16} \t{:<16} \t {:<16} \t {:<16} \t {:<18} \n".format('#0: time', '1: energy', '2: Re_L', '3: diss', '4: L', '5: U', '6: Pik'))
#output_file.write("#0: time \t 1: energy \t 2: Re_L \t 3: diss \t 4: L \t 5: U \t 6: Pik \n")

# find next block of spectra
for i in range (len(ulines)):
    line_data = ulines[i].split()
    
    try:
        spectra_block_check = line_data[1]
    except:
        pass
    
    if spectra_block_check == "Plot":
        
# calculate stats

        E = 0.0
        L = 0.0
        eps = 0.0
        Pik = 0.0
        time = float(line_data[-1])


        for k in range(1,kmax): #sum for each spectra
            line_data    = ulines[i + k].split()
            E_k          = float(line_data[1])
            if E_k < 0:
                print('Warning: Negative Energies in line {:d}'.format(i))
            T_k          = float(line_data[2])
            E           += E_k*dk[k-1]
            L           += (E_k / realk[k-1])*dk[k-1]
            eps         += (2.0 * visc * E_k * math.pow(realk[k-1],2))*dk[k-1]
            Pik         += T_k*dk[k-1]
        # final trimmings
        if E == 0:
                E = 1.0
        
        U       = math.sqrt(factor * E)
        L       *= gamma(d/2.0)*math.sqrt(math.pi)/(gamma((d+1)/2.0)*U*U) 
        Re_L    = L * U / visc

        #write data in output file

        output_file.write("{:<7.2f} \t {:<16} \t {:<16} \t{:<16} \t {:<16} \t {:<16} \t {:<18} \n".format(time,E,Re_L,eps,L,U,Pik))

#close file
output_file.close()


