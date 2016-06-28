####This will contain all the equations that will be be needed to calculate the optimal parameter space. 

import numpy as np 
import matplotlib.pyplot as plt 
from scipy import interpolate 
from scipy import integrate


####The temperatures for each region of the cryostat
Ta = 300.  #the external temperature of the cryostat (K)
Tb = 50.   #the temperature of the first stage (K)
Tc = 4.    #the temperature of the second stage (K)
Td = 1.    #the temperature of the Helium-4 stage (K)
Te = 0.35  #the temperature of the Helium-3 stage (K)
Tf = 0.1   #the temperature of the detector (servo controlled) (K)


####Common parameters
n_coax4K    = 8.             #Number of coaxial lines running to the 4-K cold plate
n_coax300mK = 4.             #Number of coaxial lines running to the 300mK stage
n_coax100mK = 4.             #Number of coaxial lines running to the 100mK stage
n_amps      = 2.             #Number of cold amplifiers installed
I_ADR       = 12.            #Maximum current fo ADR magnet (A)
sigma       = 5.67*10**(-8)  #Stefan's constant (W m^-2 K^4)
E_G10       = 28.*10**9      #Young's modulus of elasticity for G10 glass expoxy (Pa)
E_CFRP      = 120.*10**9     #Graphlite CFRP (GPa)
E_Ti15333   = 100.*10**9     #Titanium 15-3-3-3
omega_n     = 150.           #Design natural frequency for each cold mass relative to the warmer stage it is attached to (Hz)
epsilon_mli = 0.03           #the effective emissivity of the super insulation 
epsilon_al  = 0.05           #the effective emissivity of polished aluminum 
epsilon_au  = 0.02           #the effective emissivity of gold-plated surfaces
rho_cu      = 1.4            #electrical resistivity of copper 300K-70K (assumed T independant) (micro ohm cm)
rho_brass   = 4.5            #electrical resistivity of brass 70K-4K (assumed T independant) (micro ohm cm)


####Heat Conduction Integrals
#CuETB
kappa_CuETB = np.genfromtxt('thermal_tables.csv',delimiter = ',',skip_header = 1,usecols = (8,9))
kappa_CuETB_interp = interpolate.interp1d(kappa_CuETB[:,0], kappa_CuETB[:,1],fill_value = 0.,bounds_error = False)  #(W/m/K)
def K_CuETB(T):
    return integrate.quad(kappa_CuETB_interp,0,T)[0]


####Conduction along the ADR magnet leads
n_b_magnet_leads            = 4.
l_b_magnet_leads            = 660.     #(mm) ctc_ele_112
diameter_b_magnet_leads     = 0.81280  #AWG 20 wire (That number is for solid core but we actually have stranded 0.965mm diameter.)
A_b_magnet_leads            = np.pi/4.*(diameter_b_magnet_leads/1000.)**2  #(m^2)
Q_b_conduction_magnet_leads = n_b_magnet_leads*(A_b_magnet_leads/(l_b_magnet_leads/1000)*(K_CuETB(Ta) -K_CuETB(4)))
print "Q_b_conduction_magnet_leads =",Q_b_conduction_magnet_leads,"W"


####Joule heating in the ADR magnet leads --> This needs to be updated with just the copper part on the 4K stage. I assume that joule heating is not important on the superconducting leads
R_b_single_lead    = (rho_cu/100.)*(l_b_magnet_leads/1000.)/A_b_magnet_leads  #resistance of single lead (micro Ohms)
R_b_single_contact = 100.                                                     #(micro Ohms)
R_b_magnet_leads   = n_b_magnet_leads*(R_b_single_lead+R_b_single_contact)    #(micro Ohms)
Q_b_joule_heating  = (I_ADR/2)**2*R_b_magnet_leads/(10**6)                    #divided current by two since there are two wires for each lead
print "Q_b_joule_heating =",Q_b_joule_heating,"W"
