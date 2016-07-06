#%matplotlib inline
import argparse
import csv
import itertools
import numpy as np
import matplotlib.pyplot as plt


''' This code will calculate the optimal parameter space for our Helmholtz coil.
 It will automatically generate a csv file when the current through each solenoid is specified in the command-line. '''


####The temperatures for each region of the cryostat
Ta = 300.  #The external temperature of the cryostat (K)
Tb = 50.   #The temperature of the first stage (K)
Tc = 4.    #The temperature of the second stage (K)
Td = 1.    #The temperature of the Helium-4 stage (K)
Te = 0.35  #The temperature of the Helium-3 stage (K)
Tf = 0.1   #The temperature of the detector (servo controlled) (K)


####Common parameters
#rho_cu = 1.4 * 10**-6     #Electrical resistivity of copper 300K-70K (assumed T independent) (Ohm cm) --> Where did they get this number?
rho_cu = 1.68 * 10**-8     #Electrical resistivity of copper (Ohm m))
mu_0 = 4 * np.pi * 10**-7  #Permeability constant (T m / A)
mu_r = 1.                  #Relative permeability constant (T m / A) ---> This will need to be experimentally determined for our cryostat.


parser = argparse.ArgumentParser(description='Type in a value for current (A).')
parser.print_help()
parser.add_argument('I', type=float, help='Current (A)')
args = parser.parse_args()


####Important parameters
r = 0.0254 * 13.                                     #Radius of each solenoid (m)
r_in = r / 0.0254									 #Radius of each solenoid (in)
B_0 = 65. * 10**-5.                                  #10x magnetic field of Earth (T)
I = args.I                                           #Current through each solenoid (A)
N = B_0 * 5 * np.sqrt(5) * r / (8 * mu_0 * mu_r * I) #Number of turns in each solenoid
l = 2 * np.pi * r * N                                #Length of wire for each solenoid (m)
l_in = l / 0.0254 / 12.								 #Length of wire for each solenoid (ft)

print "\nr =", r_in, "in"
print "l =", l_in, "ft"
print "B_0 =", B_0, "T"
print "I =", I, "A"
print "N =", N


awg_all = np.genfromtxt('magnetwiretable.csv',delimiter = ',',skip_header = 1,usecols = (0,1,3,4))  #This will import the data for gauge, cross-sectional area, cost, and weight in that order.


def R(l, A):
	''' Calculates the resistance of each solenoid of the Helmholtz coil. '''
	return rho_cu * l / (A * 10**-6)  #(Ohm)


def L(N, A, l):
	''' Calculates the inductance of each solenoid of the Helmholtz coil. '''
	return mu_0 * N**2 * A * 10**-6 / l  #(H)

  
def Q_joule_heating(I, R):
	''' Calculates the Joule heating for our Helmholtz coil. '''
	return I**2 * R  #(W)

	
with open('helmholtz4.csv', 'wb') as csvfile:
	writer = csv.writer(csvfile)
	writer.writerow(['AWG', 'Cross-sectional Area (mm^2)', 'Cost (per ft)', 'Weight (lb per ft)', 'Resistance (Ohm)', 'Inductance (H)', 'Magnetic Field (T)', 'Joule Heating (W)'])
	for a, b, c, d in itertools.izip(awg_all[:,0], awg_all[:,1], awg_all[:,2], awg_all[:,3]):
		awg = a
		A = b
		cost = c * l_in
		weight = d * l_in
		resistance = R(l, A)
		inductance = L(N, A, l)
		magnetf = B_0
		heating = Q_joule_heating(I, resistance)
		writer.writerow([awg, A, cost, weight, resistance, inductance, magnetf, heating])

csvfile.close()

newdata = np.genfromtxt('helmholtz4.csv',delimiter = ',',skip_header = 1,usecols = (0,3,7))  #This will import the data from the csv file just created.

awg = newdata[:,0]
weight = newdata[:,1]
heating = newdata[:,2]

#Joule Heating vs. Gauge
plot1 = plt.plot(awg, heating) 
plt.xlabel('Gauge (AWG)')
plt.ylabel('Joule Heating (W)')
plt.title("Joule Heating vs. Gauge")
plt.show()

#Weight vs. Gauge
plot2 = plt.plot(awg, weight) 
plt.xlabel('Gauge (AWG)')
plt.ylabel('Weight (lb)')
plt.title("Weight vs. Gauge")
plt.show()

# Nice Table of Inputs to notebook and PDF

'''print '_____________________________________________________________________________________'
print 'TABLE OF INPUTS AND RESULTANT DETECTOR PARAMETERS'
print '_____________________________________________________________________________________'
print 'Parameter     Value     Units             Type         Description'
print '_____________________________________________________________________________________'
print 'f_o           %5.2e  Hz                Measured     Resonant frequency' % (omega/2.0/np.pi)
print 'T_c           %5.2e  K                 Measured     Critical temperatures' % T_c
print 'Qi            %5.2e  Dimensionless     Measured     Internal Q' % Q_i
print 'Qc            %5.2e  Dimensionless     Measured     Coupling Q' % Q_c
print 'Qr            %5.2e  Dimensionless     Measured     Total Q' % Q_r
print 'eta_a         %5.2e  Dimensionless     Assumed      Microwave coupling efficiency' % eta_a
print 'eta_o         %5.2e  Dimensionless     Assumed      Optical coupling efficiency' % eta_o
print 'tau_max       %5.2e  s                 Guess        Maximum thermal time constant' % tau_max
print 'tau_qp        %5.2e  s                 Measured     QP time constant' % tau_qp
print 'V             %5.2e  um^3              Estimated    Detector active volume' % V
print '2*Delta_0     %5.2e  eV                Input        QP pair binding energy' % (2*Delta_0/1.602e-19)
print 'n_star        %5.2e  um^-3             Assumed      "Crossover" density of states' % n_star
print 'beta          %5.2e  Dimensionless     Calculated   Ratio freq/amp response' % (bbeta)
print 'chi_c         %5.2e  Dimensionless     Assumed      Coupling Q' % (chi_c)
print 'chi_qp        %5.2e  Dimensionless     Assumed      Frac. resonator diss. from QP' % (chi_qp)
print 'S_TLS         %5.2e  Hz^-1             Assumed      TLS amplitude' % (S_TLS)
print 'T_a           %5.2e  K                 Assumed      Amplifier noise temperature' % T_a
print 'T_BB          %5.2e  K                 Assumed      Blackbody temperature' % (T_BB)
print 'n_0           %5.2e  Dimensionless     Calculated   Photon occupation number' % (n_0)
print 'Omega_lens    %5.2e  steradians        Calculated   Solid angle of blackbody' % (Omega_lens)
print 'A_lenslet     %5.2e  m^2               Assumed      Effective area of lens' % (A_lenslet)
print 'P_o           %5.2e  W                 Calculated   Blackbody power on detector' % (P_o)
print 'Photon rate   %5.2e  s^-1              Approximate  Photon arrival rate @ detector' % (P_o/h/nu)
print 'P_a           %5.2e  W                 P_a = P_o    Microwave power (assumed)' % (P_a)
print 'Gamma         %5.2e  s^-1              Calculated   Total QP generation rate' % (GGamma)
print 'QP prs/photon %5.2e  Dimensionless     Gamma/2      # QP pairs per submm photon' % ((GGamma/2.0)/(P_o/h/nu))
print 'n_qp          %5.2e  um^-3             Calculated   QP number density w/o thermal' % (n_qp)
print 'n_qp_therm    %5.2e  um^-3             Calculated   Therm QP number density' % (n_qp_therm/1.0e18)
print 'N_qp          %5.2e  Dimensionless     Calculated   Number of QPs' % (N_qp)
print '______________________________________________________________________________________'

outfile = open('KID_Parameters.txt', 'w')
#outfile.write('%5.1f' % column)
outfile.write('_____________________________________________________________________________________\n')
outfile.write('TABLE OF INPUTS AND RESULTANT DETECTOR PARAMETERS\n')
outfile.write('_____________________________________________________________________________________\n')
outfile.write('Parameter     Value     Units             Type         Description\n')
outfile.write('_____________________________________________________________________________________\n')
outfile.write('f_o           %5.2e  Hz                Measured     Resonant frequency\n' % (omega/2.0/np.pi))
outfile.write('T_c           %5.2e  K                 Measured     Critical temperature\n' % T_c)
outfile.write('Qi            %5.2e  Dimensionless     Measured     Internal Q\n' % Q_i)
outfile.write('Qc            %5.2e  Dimensionless     Measured     Coupling Q\n' % Q_c)
outfile.write('Qr            %5.2e  Dimensionless     Measured     Total Q\n' % Q_c)
outfile.write('eta_a         %5.2e  Dimensionless     Assumed      Microwave coupling efficiency\n' % eta_a)
outfile.write('eta_o         %5.2e  Dimensionless     Assumed      Optical coupling efficiency\n' % eta_o)
outfile.write('tau_max       %5.2e  s                 Guess        Maximum thermal time constant\n' % tau_max)
outfile.write('tau_qp        %5.2e  s                 Measeured    QP time constant\n' % tau_qp)
outfile.write('V             %5.2e  um^3              Assumed      Detector active volume\n' % V)
outfile.write('2*Delta_0     %5.2e  eV                Input        QP pair binding energy\n' % (2*Delta_0/1.602e-19))
outfile.write('n_star        %5.2e  um^-3             Assumed      "Crossover" density of states\n' % n_star)
outfile.write('beta          %5.2e  Dimensionless     Calculated   Ratio freq/amp response\n' % (bbeta))
outfile.write('chi_c         %5.2e  Dimensionless     Assumed      Coupling Q\n' % (chi_c))
outfile.write('chi_qp        %5.2e  Dimensionless     Assumed      Frac. resonator diss. from QP\n' % (chi_qp))
outfile.write('S_TLS         %5.2e  Hz^-1             Assumed      TLS amplitude\n' % (S_TLS))
outfile.write('T_a           %5.2e  K                 Assumed      Amplifier noise temperature\n' % T_a)
outfile.write('T_BB          %5.2e  K                 Assumed      Blackbody temperature\n' % (T_BB))
outfile.write('n_0           %5.2e  Dimensionless     Calculated   Photon occupation number\n' % (n_0))
outfile.write('Omega_lens    %5.2e  steradians        Calculated   Solid angle of blackbody\n' % (Omega_lens))
outfile.write('A_lenslet     %5.2e  m^2               Assumed      Effective area of lens\n' % (A_lenslet))
outfile.write('P_o           %5.2e  W                 Calculated   Blackbody power on detector\n' % (P_o))
outfile.write('Photon rate   %5.2e  s^-1              Approximate  Photon arrival rate @ detector\n' % (P_o/h/nu))
outfile.write('P_a           %5.2e  W                 P_a = P_o    Microwave power (assumed)\n' % (P_a))
outfile.write('Gamma         %5.2e  s^-1              Calculated   Total QP generation rate\n' % (GGamma))
outfile.write('QP prs/photon %5.2e  Dimensionless     Gamma/2      # QP pairs per submm photon\n' % ((GGamma/2.0)/(P_o/h/nu)))
outfile.write('n_qp          %5.2e  um^-3             Calculated   QP number density w/o thermal\n' % (n_qp))
outfile.write('n_qp_therm    %5.2e  um^-3             Calculated   Therm QP number density' % (n_qp_therm/1.0e18))
outfile.write('N_qp          %5.2e  Dimensionless     Calculated   Number of QPs\n' % (N_qp))
outfile.write('_____________________________________________________________________________________\n')'''

print "\nOutput has been written to helmholtz4.csv."
