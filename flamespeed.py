"""
Sensitivity analysis for a freely-propagating, premixed methane-air
flame. Computes the sensitivity of the laminar flame speed with respect
to each reaction rate constant.
"""



import sys
import csv
import cantera as ct
import numpy as np
from matplotlib.pylab import *

#Mechanisms used for the process
# IdealGasMix object used to compute mixture properties
gas = ct.Solution('gri30.xml', 'gri30_mix')

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'CH4:0.45, O2:1.0, N2:3.76'
npoints = 1
initial_grid = np.linspace(0, 0.03, 5)  # m

# Flame speed storage variables
Sit = np.zeros(npoints,'d')
Sip = np.zeros(npoints, 'd')
Ti = np.zeros(npoints, 'd')
Pi = np.zeros(npoints, 'd')
print()
print('\n Sl(T), p = 1 atm, fi = 1\n')

for k in range(npoints):
    gas.TPX = Tin, p, reactants
    f = ct.FreeFlame(gas, initial_grid)

    # Solve with the energy equation disabled
    f.energy_enabled = False
    f.set_max_jac_age(10, 10)
    f.set_time_step(1e-5, [2, 5, 10, 20])
    f.solve(loglevel=1, refine_grid=False)

    # Solve with the energy equation enabled
    f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
    f.energy_enabled = True
    f.solve(loglevel=1, refine_grid=True)

    Sit[k] = f.u[0]
    Ti[k]= Tin
    print('\nmixture-averaged flamespeed versus temp= {:7f} m/s\n'.format(f.u[0]))

    Tin += 50


for k in range(npoints):
    Tin = 300
    gas.TPX = Tin, p, reactants
    f = ct.FreeFlame(gas, initial_grid)

    # Solve with the energy equation disabled
    f.energy_enabled = False
    f.set_max_jac_age(10, 10)
    f.set_time_step(1e-5, [2, 5, 10, 20])
    f.solve(loglevel=1, refine_grid=False)

    # Solve with the energy equation enabled
    f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
    f.energy_enabled = True
    f.solve(loglevel=1, refine_grid=True)

    Sip[k] = f.u[0]
    Pi[k]= p
    print('\nmixture-averaged flamespeed versus pressure = {:7f} m/s\n'.format(f.u[0]))

    p+=50000
'''
for i in range(npoints):
    print ("\nSi = {:7f} m/s\n".format(Si[i]))
    print ("For "+str(Ti[i])+' K laminar flame speed = '+str(Si[i])+" m\s\n")


# create plot T
plot(Ti,Sit, '-', color = 'pink')
xlabel(r'Temp [K]', fontsize=20)
ylabel("Sl [m\s]")
title(r'Laminar flame speed of $H_{2}$ + Air mixture at $\Phi$ = 1, and P = 1 atm',
fontsize=22,horizontalalignment='center')
axis([290,600,0.32,1.5])
grid()
#show()
savefig('Laminar_T.png', bbox_inches='tight')

#create plot p
plot(Pi,Sip, '-', color = 'pink')
xlabel(r'p [Pa]', fontsize=20)
ylabel("Sl [m\s]")
title(r'Laminar flame speed of $H_{2}$ + Air mixture at $\Phi$ = 1, and P = 1 atm',
fontsize=22,horizontalalignment='center')
axis([101325,900000,0.32,1.5])
grid()
#show()
savefig('Laminar_p.png', bbox_inches='tight')
'''