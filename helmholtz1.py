import csv
import collections
import numpy as np 
#import matplotlib.pyplot as plt


''' This code contains all the equations needed to calculate the optimal parameter space for our Helmholtz coil.
	 It will automatically generate a csv file of values when the length (l), number of turns in each solenoid (N),
	 current through each solenoid (I), and radius of each solenoid (r) are specified in the command-line (?). '''


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

####Command-line parameters
#Give me l, N, I, r
l = 1.	 #Length of wire (m)
N = 20.  #Number of turns in each solenoid
I = 1.	 #Current through each solenoid (A)
r = 0.5  #Radius of each solenoid (m)


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
	''' Calculates the inductance of each solenoid of the Helmholtz coil. '''  #Should we bother finding the relaxation time constant tau?
	return mu_0 * N**2 * A * 10**-6 / l  #(H)


def B_0(N, I, r):
	''' Calculates the magnetic field generated by a Helmholtz coil at its center. '''
	#return (np.sqrt((4/5)**3) * mu_0 * mu_r * N * I) / r  #(T)
	return 8 * mu_0 * mu_r * N * I / (5 * np.sqrt(5) * r)  #(T)
	
  
def Q_joule_heating(I, R):
	''' Calculates the Joule heating for our Helmholtz coil. '''
	return I**2 * R  #(W)


with open('helmholtz4.csv', 'wb') as csvfile:
	writer = csv.writer(csvfile)
	writer.writerow(['AWG', 'Cross-sectional Area (mm^2)', 'Resistance (Ohm)', 'Inductance (H)', 'Magnetic Field (T)', 'Joule Heating (W)'])
	for key, value in awg_od.iteritems():
		A = value
		resistance = R(l, A)
		inductance = L(N, A, l)
		magnetf = B_0(N, I, r)
		heating = Q_joule_heating(I, resistance)
		writer.writerow([key, A, resistance, inductance, magnetf, heating])

csvfile.close()

# I don't know what to do with this --> V = L * dI / dt

''' ####Joule heating in the Helmholtz coil
# Joule heating (H) = I**2 * R (current^2 * resistance of copper) in W
Q_b_joule_heating  = (I_ADR/2)**2*R_b_magnet_leads/(10**6)  #divided current by two since there are two wires for each lead
print "Q_b_joule_heating =",Q_b_joule_heating,"W" '''
