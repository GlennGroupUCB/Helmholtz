import argparse
import csv
import collections
import numpy as np 
from termcolor import colored, cprint
#import matplotlib.pyplot as plt


''' This code will calculate the optimal parameter space for our Helmholtz coil.
	 It will automatically generate a csv file when the length of wire (l), number of turns in each solenoid (N),
	 current through each solenoid (I), and radius of each solenoid (r) are specified in the command-line. '''


####The temperatures for each region of the cryostat
Ta = 300.  #The external temperature of the cryostat (K)
Tb = 50.   #The temperature of the first stage (K)
Tc = 4.	   #The temperature of the second stage (K)
Td = 1.	   #The temperature of the Helium-4 stage (K)
Te = 0.35  #The temperature of the Helium-3 stage (K)
Tf = 0.1   #The temperature of the detector (servo controlled) (K)


####Common parameters
#rho_cu = 1.4 * 10**-6	   #Electrical resistivity of copper 300K-70K (assumed T independent) (Ohm cm) --> Where did they get this number?
rho_cu = 1.68 * 10**-8	   #Electrical resistivity of copper (Ohm m))
mu_0 = 4 * np.pi * 10**-7  #Permeability constant (T m / A)
mu_r = 1.				   #Relative permeability constant (T m / A) ---> This will need to be experimentally determined for our cryostat.


#parser = argparse.ArgumentParser(description='''Type in values for length of wire (m), number of turns, current (A), 
#and radius (m) in that order.''')
''' parser.print_help()
parser.add_argument('l', type=float, help='Length of wire (m)')
parser.add_argument('N', type=float, help='Number of turns')
parser.add_argument('I', type=float, help='Current (A)')
parser.add_argument('r', type=float, help='Radius (m)')
args = parser.parse_args() '''


''' l = args.l
N = args.N
I = args.I
r = args.r '''


####Important parameters
r = 0.0254 * 13.                                     #Radius of each solenoid (m)
l = 2 * np.pi * r	                                 #Length of wire (m)
B_0 = 65. * 10**-5.                                  #10x magnetic field of Earth (T)
I = 1.	                                             #Current through each solenoid (A)
N = B_0 * 5 * np.sqrt(5) * r / (8 * mu_0 * mu_r * I) #Number of turns in each solenoid

print "r =", r / 0.0254, "in"
print "l =", l / 0.0254 / 12., "ft"
print "B_0 =", B_0, "T"
print "I =", I, "A"
print "N =", N


awg = {'4_0': 107., '3_0': 85., '2_0': 67.4, '1_0': 53.5, '1': 42.4, '2': 33.6, '3': 26.7,
		'4': 21.2, '5': 16.8, '6': 13.3, '7': 10.5, '8': 8.37, '9': 6.63, '10': 5.26, '11': 4.17,
		'12': 3.31, '13': 2.62, '14': 2.08, '15': 1.65, '16': 1.31, '17': 1.04, '18': 0.823,
		'19': 0.653, '20': 0.518, '21': 0.410, '22': 0.326, '23': 0.258, '24': 0.205, '25': 0.162,
		'26': 0.129, '27': 0.102, '28': 0.0810, '29': 0.0642, '30': 0.0509, '31': 0.0404, '32': 0.0320,
		'33': 0.0254, '34': 0.0201, '35': 0.0160, '36': 0.0127, '37': 0.0100, '38': 0.00797,
		'39': 0.00632, '40': 0.00501}  #American wire gauge (AWG) cross sectional area (mm^2)

awg_od = collections.OrderedDict(sorted(awg.items(), key=lambda t: t[1]))  #This orders awg in descending order according to the key.



def R(l, A):
	''' Calculates the resistance of each solenoid of the Helmholtz coil. '''
	return rho_cu * l / (A * 10**-6)  #(Ohm)


def L(N, A, l):
	''' Calculates the inductance of each solenoid of the Helmholtz coil. '''
	return mu_0 * N**2 * A * 10**-6 / l  #(H)

  
def Q_joule_heating(I, R):
	''' Calculates the Joule heating for our Helmholtz coil. '''
	return I**2 * R  #(W)


with open('helmholtz3.csv', 'wb') as csvfile:
	writer = csv.writer(csvfile)
	writer.writerow(['AWG', 'Cross-sectional Area (mm^2)', 'Resistance (Ohm)', 'Inductance (H)', 'Magnetic Field (T)', 'Joule Heating (W)'])
	for key, value in awg_od.iteritems():
		A = value
		resistance = R(l, A)
		inductance = L(N, A, l)
		magnetf = B_0
		heating = Q_joule_heating(I, resistance)
		writer.writerow([key, A, resistance, inductance, magnetf, heating])

csvfile.close()

print colored("\nOutput has been written to helmholtz3.csv.", 'green')
