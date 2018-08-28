# -*- coding: utf-8 -*-
# Discrete Monte Carlo Model for Secondary Electron Emission (SEE)
# Hsing-Yin (Iris) Chang (irischang@ucla.edu)
# Dpt. of Materials Science & Engineering (MSE) at UCLA
# Jaime Marian Group
# current edited date Aug 27, 2018

from __future__ import division
import sys
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')
from numpy import *
import numpy as np
import scipy.interpolate
from scipy.interpolate import PchipInterpolator
import random
from datetime import datetime


if len(sys.argv) != 4:
	print "The Single Scattering Monte Carlo Model"
	exit("Usage: %s [primary electron energy (keV)] [incident angle (radian)] [sampling times] " % sys.argv[0])

input1 = float(sys.argv[1])                  # primary electron energy [keV]
input2 = float(sys.argv[2])                  # incident angle [radian]
input3 = int(sys.argv[3])                    # sampling times


########## Hexagonal Boron Nitride (h-BN) ##########
# source: wikipedia & http://www.ioffe.ru/SVA/NSM/Semicond/BN/basic.html#hexagonal
Z = 12.0                                     # Z = Z_B+Z_N
Na = 6.022E23                                # Avagadro's number
rho = 2.18                                   # bulk density [g/cm^3]
AW = 24.8177                                 # molecular weight [g/mol]
N = 5.29E22                                  # volume number density [molecules/cm^3] = rho*Na/AW
EA = 4.5                                     # electron affinity [eV]
Eg = 5.2                                     # bandgap energy [eV]
Nv = 4.0                                     # number of valence electrons per atom (for elemental solids) or molecule (for compounds)
Ep = 28.8*np.sqrt(Nv*rho/AW)                 # free-electron plasmon frequency [eV]

########## Coefficient ########## 
beta = -0.1+0.944*(Ep**2+Eg**2)**(-0.5)+0.069*rho**0.1
gamma = 0.191*rho**(-0.5)
U = Nv*rho/AW                                # (= Ep**2/829.4)
C = 1.97-0.91*U
D = 53.4-20.8*U


########## Physical constants ########## 
e = 4.8032068E-10              # electron charge [esu]
m_e = 9.10938E-28              # electron mass [g]
c = 2.99792458E10              # speed of light [cm/s]
r_e = 2.8179403E-13            # electron radius [cm]
hbar = 1.0546E-27              # reduced Planck constant [erg*s] or [cm^2-g/s]
a_0 = 0.529177                 # Bohr radius [Angstrom]


# A binary tree node
class Node:
# Constructor to create a new node
	def __init__(self, E, x, y, z, theta, psi, steplength):
		self.E = E
		self.x = x
		self.y = y
		self.z = z
		self.theta = theta
		self.psi = psi
		self.steplength = steplength
		self.left = None
		self.right = None


def buildTree(root, height, count):
	data = np.asarray([1000*(root.E), root.x, root.y, root.z]).reshape(1,4)
	if height == 0:
		# print "Terminated."
		return count
	elif root.E <= EA*1E-3:
		# print "The energy of secondry electron falls below surface potential barrier."
		return count
	elif root.z < 0:
		# print "Secondry electron emitted out of the metal surface."
		count_out = count+1
		count = count_out
		data = np.asarray([1000*(root.E)-EA, root.x, root.z, root.theta]).reshape(1,4)
		print data
		return count

	mfp_el = mfp_elastic(root.E)
	# print "elsetic mean free path: ", mfp_el
	mfp_in = mfp_inelastic(root.E)
	# print "elsetic mean free path: ", mfp_el
	mfp_ph = mfp_phonon(root.E)
	# print 'phonon mean free path', mfp_ph
	mfp_pol = mfp_polaron(root.E)
	# print 'polaron mean free path', mfp_pol
	mfp_T = 1/((1/mfp_el)+(1/mfp_in)+(1/mfp_ph)+(1/mfp_pol))
	# print "total mean free path: ", mfp_T
	Test_scatter = np.random.uniform(0,1) 
	P_el = mfp_T/mfp_el
	# print 'P_el', P_el
	if Test_scatter > P_el:
		#print 'inelastic scattering'
		R_E = np.random.uniform(0,1)                    # random number
		if Test_scatter < (1/mfp_el+1/mfp_in)/(1/mfp_T):
			#print 'inelastic'
			left = inelasticScatter(root, root.E, root.x, root.y, root.z, root.theta, root.psi, root.steplength, True, R_E)
			right = inelasticScatter(root, root.E, root.x, root.y, root.z, root.theta, root.psi, root.steplength, False, R_E)
			root.left = left
			root.right = right
			count_left = buildTree(left, height-1, count)
			count = count_left
			count_right = buildTree(right, height-1, count)
			count_leftandright = count_right
			count = count_leftandright
			return count
		elif Test_scatter < (1/mfp_el+1/mfp_in+1/mfp_ph)/(1/mfp_T):
			#print 'phonon'
			next = phononExcitation(root, root.E, root.x, root.y, root.z, root.theta, root.psi, root.steplength)
			root.left = next
			root.right = None
			count_left = buildTree(next, height-1, count)
			count = count_left
			return count
		else:
			#print 'polaron'
			return count
	else:
		#print 'elastic'
		next = elasticScatter(root, root.E, root.x, root.y, root.z, root.theta, root.psi, root.steplength)
		root.left = next
		root.right = None
		count_left = buildTree(next, height-1, count)
		count = count_left
		return count


A1 = []
B1 = []
b1 = []
x3 = []

# Define eV values at which we have extracted data for
x1 = np.linspace(50,100,6)
x2 = np.linspace(150,1000,18)
x3 = np.concatenate((x1, x2), axis=0)

# Define eV values that we wish to interpolate
X = np.linspace(50,1000, 95100)

arr_A1 = np.empty([24,181])
# Extract data from .txt file
with open("../sample_inputs/DESCS_hBN.txt", "r+") as DESCS: 
	for line in DESCS:
		data1 = line.split()
		A1.append(data1)
arr_A1 = np.asarray(A1).astype(np.float)

arr_B1 = np.empty([95100,181])
# For each column, interpolate
for i in range(181):
	bi1 = PchipInterpolator(x3, arr_A1[:,i])
	b1 = bi1(X)
	B1.append(b1)
arr_B1 = np.asarray(np.transpose(B1)).astype(np.float)


def mfp_elastic(E):
	m = 176.0
	q = 0.733
	alpha_c = 0.5
	R_c = np.tanh(alpha_c*(1000*E/Eg)**2)
	sigma_el = R_c*m/((1000*E)**q)*1E-16
	return 2*AW/(Na*rho*sigma_el)*1E8


########## ELF ##########
# x - energy loss, y = momentum transfer 
def ELF_hBN(x,y):
	alpha = 0.05
	def E1(y):
		return 8.65+alpha*hbar**2*(y*np.sqrt(1.6*10**(-19))/hbar)**2/(2*m_e)*6.242*10**18
	def E2(y):
		return 19.0+alpha*hbar**2*(y*np.sqrt(1.6*10**(-19))/hbar)**2/(2*m_e)*6.242*10**18
	def E3(y):
		return 25.2+alpha*hbar**2*(y*np.sqrt(1.6*10**(-19))/hbar)**2/(2*m_e)*6.242*10**18
	def E4(y):
		return 35.8+alpha*hbar**2*(y*np.sqrt(1.6*10**(-19))/hbar)**2/(2*m_e)*6.242*10**18
	def E5(y):
		return 51.0+alpha*hbar**2*(y*np.sqrt(1.6*10**(-19))/hbar)**2/(2*m_e)*6.242*10**18
	def E6(y):
		return 65.0+alpha*hbar**2*(y*np.sqrt(1.6*10**(-19))/hbar)**2/(2*m_e)*6.242*10**18
	gamma1 = 0.5
	gamma2 = 5.0
	gamma3 = 9.5
	gamma4 = 10.0
	gamma5 = 20.0
	gamma6 = 20.0
	A1 = 6.6
	A2 = 18.0
	A3 = 270.0
	A4 = 140.0
	A5 = 80.0
	A6 = 20.0
	return A1*gamma1*x/((E1(y)**2-x**2)**2+(gamma1*x)**2)+A2*gamma2*x/((E2(y)**2-x**2)**2+(gamma2*x)**2)+A3*gamma3*x/((E3(y)**2-x**2)**2+(gamma3*x)**2)+A4*gamma4*x/((E4(y)**2-x**2)**2+(gamma4*x)**2)+A5*gamma5*x/((E5(y)**2-x**2)**2+(gamma5*x)**2)+A6*gamma6*x/((E6(y)**2-x**2)**2+(gamma6*x)**2)


########## Integrating IMFP ########## 
def mfp_inelastic(E):
	list_B = []
	for i in range(1000):
	# Iterate inelastic mean free path equation for all Ep_cgs, n, kv
	# Define equation to integrate
		Ep_cgs = 1000*E
		x = np.linspace(1,1000,1000)
		A = (m_e*e**2)/(2*np.pi*hbar**2*Ep_cgs)   # Constant term
		B = (ELF_hBN(i,0))*((1.0-(x[i]/Ep_cgs))*np.log(4.0/(x[i]/Ep_cgs))-(7.0/4.0)*(x[i]/Ep_cgs)+(x[i]/Ep_cgs)**(3.0/2.0)-(33.0/32.0)*(x[i]/Ep_cgs)**2)	
		list_B.append(B)
	array_B = np.asarray(list_B).astype(np.float)
	array_x = x[int(Eg):int(1000*E/2)]
	array_B = array_B[int(Eg):int(1000*E/2)]      # Limit integrating function to integral bounds
	s = np.trapz(array_B, x=[array_x])+1E-15
	return (1.0/(A*s))*1E8


def mfp_phonon(E):
	r_Bohr = 0.529                                # Bohr radius [A]
	epsilon_0 = 5.09
	epsilon_inf = 4.1
	dE = 0.1                                      # electron energy loss [eV]
	k_B = 8.61733035E-5                           # Boltzmann constant [eV/K]
	T = 300.0                                     # occupation number for the phonon level at temperature T
	n_occ = 1.0/(np.exp(dE/(k_B*T))-1.0)
	mfp_phonon_inv = 1.0/r_Bohr*((epsilon_0-epsilon_inf)/(epsilon_0*epsilon_inf))*(dE/(1000*E))*0.5*(n_occ+1.0)*np.log((1.0+np.sqrt(1.0-dE/(1000*E)))/(1.0-np.sqrt(1.0-dE/(1000*E))))
	return 1.0/mfp_phonon_inv                     # phonon mean free path [A]


def mfp_polaron(E):
	C = 1.0                                       # [nm^-1]
	gamma = 0.1                                   # [eV^-1]
	return 1.0/(C*np.exp(-gamma*(1000*E)))*10     # polaron mean free path [A]

A_ThetaEl = []
B_ThetaEl = []
b_ThetaEl = []
x3_ThetaEl = []

# Define eV values at which we have extracted data for
x1_ThetaEl = np.linspace(50,100,6)
x2_ThetaEl = np.linspace(150,1000,18)
x3_ThetaEl = np.concatenate((x1_ThetaEl, x2_ThetaEl), axis=0)

# Define eV values that we wish to interpolate
X_ThetaEl = np.linspace(50,1000,95001)

arr_A_ThetaEl = np.empty([24,100])
# Extract data from .txt file
with open("../sample_inputs/ThetaEl_hBN.txt", "r+") as ESA: 
	for line in ESA:
		data = line.split()
		A_ThetaEl.append(data)
arr_A_ThetaEl = np.asarray(A_ThetaEl).astype(np.float) # Convert to float array
arr_A_ThetaEl = np.transpose(arr_A_ThetaEl)

# For each column, interpolate
for i in range(100):
	bi_ThetaEl = PchipInterpolator(x3_ThetaEl, arr_A_ThetaEl[:,i])
	b_ThetaEl = bi_ThetaEl(X_ThetaEl)
	B_ThetaEl.append(b_ThetaEl)
arr_B_ThetaEl = np.asarray(B_ThetaEl).astype(np.float)


def elasticScatter(root, E, x, y, z, theta, psi, steplength): 
	mfp_el = mfp_elastic(E)
	mfp_in = mfp_inelastic(E)
	mfp_ph = mfp_phonon(E)
	mfp_pol = mfp_polaron(E)
	mfp_T = 1/((1/mfp_el)+(1/mfp_in)+(1/mfp_ph)+(1/mfp_pol))
	R1 = np.random.uniform(0,1)
	steplength = -mfp_T*np.log(R1)
	idx = int(100*round(1000*E,2)-4999)
	R2 = np.random.uniform(0,1)
	theta_PE = np.deg2rad(arr_B_ThetaEl[int(100*round(R2,2)-1)][idx])
	R3 = np.random.uniform(0,1)
	psi_PE = 2*np.pi*R3
	theta_new = np.arccos(np.cos(theta)*np.cos(theta_PE)+np.sin(theta)*np.sin(theta_PE)*np.cos(psi_PE))
	psi_new = np.arcsin(np.sin(theta_PE)*np.sin(psi_PE)/np.sin(theta_new))+psi
	x = x+steplength*(np.sin(theta_new)*np.cos(psi_new))
	y = y+steplength*(np.sin(theta_new)*np.sin(psi_new))
	z = z+steplength*(np.cos(theta_new))
	theta = theta_new
	psi = psi_new
	return Node(E, x, y, z, theta, psi, steplength)


def phononExcitation(root, E, x, y, z, theta, psi, steplength):
	E_prime = E-(0.1/1000)
	B = (1000*E+1000*E_prime+2*np.sqrt(1000*E*1000*E_prime))/(1000*E+1000*E_prime-2*np.sqrt(1000*E*1000*E_prime))
	R1 = np.random.uniform(0,1)           # random number	
	theta_PE = np.arccos(((1000*E+1000*E_prime)/(2*np.sqrt(1000*E*1000*E_prime)))*(1-B**R1)+B**R1)
	R2 = np.random.uniform(0,1)
	psi_PE = 2*np.pi*R2
	mfp_el = mfp_elastic(E)
	mfp_in = mfp_inelastic(E)
	mfp_ph = mfp_phonon(E)
	mfp_pol = mfp_polaron(E)
	mfp_T = 1/((1/mfp_el)+(1/mfp_in)+(1/mfp_ph)+(1/mfp_pol))
	R_r = np.random.uniform(0,1)
	steplength = -mfp_T*np.log(R_r)
	theta_new = np.arccos(np.cos(theta)*np.cos(theta_PE)+np.sin(theta)*np.sin(theta_PE)*np.cos(psi_PE))
	psi_new = np.arcsin(np.sin(theta_PE)*np.sin(psi_PE)/np.sin(theta_new))+psi
	x = x+steplength*(np.sin(theta_new)*np.cos(psi_new))
	y = y+steplength*(np.sin(theta_new)*np.sin(psi_new))
	z = z+steplength*(np.cos(theta_new))
	theta = theta_new
	psi = psi_new
	E = E_prime
	return Node(E, x, y, z, theta, psi, steplength)


list_A = []
list_B = []
b = []
t3 = []


# Define eV values at which we have extracted data for
t1 = np.linspace(5.2,5.2,1)
t2 = np.linspace(50,1000,24)
t3 = np.concatenate((t1, t2), axis=0)
# Define eV values that we wish to interpolate
t = np.linspace(5.2,1000,99481)

# Extract data from .txt file
with open("../sample_inputs/ELoss_hBN.txt", "r+") as ELoss: 
	for line in ELoss:
		data = line.split()
		list_A.append(data)
arr_A = np.asarray(list_A).astype(np.float) #Convert to float array
arr_A = np.transpose(arr_A)

# For each column, interpolate
for i in range(100):
	bi = PchipInterpolator(t3, arr_A[:,i])
	b = bi(t)
	list_B.append(b)
arr_B = np.asarray(list_B).astype(np.float)
arr_C = np.asarray([[0]*519]*100).astype(np.float)
arr_D = np.concatenate((arr_C,arr_B), axis=1)


def inelasticScatter(root, E, x, y, z, theta, psi, steplength, isGenerateLeft, R_E):
	R_s = random.uniform(0,1)
	idx = int(100*round(1000*E,2)-1)
	W = arr_D[int(100*round(R_s,2)-1)][idx]*10**(-3)
	theta_PE = np.arcsin(np.sqrt(W/E))
	R3 = np.random.uniform(0,1)
	psi_PE = 2*np.pi*R3
	if isGenerateLeft:
		E = E-W
		R_r = np.random.uniform(0,1)
		mfp_el = mfp_elastic(E)
		mfp_in = mfp_inelastic(E)
		mfp_ph = mfp_phonon(E)
		mfp_pol = mfp_polaron(E)
		mfp_T = 1/((1/mfp_el)+(1/mfp_in)+(1/mfp_ph)+(1/mfp_pol))
		steplength = -mfp_T*np.log(R_r)
		theta_new = np.arccos(np.cos(theta)*np.cos(theta_PE)+np.sin(theta)*np.sin(theta_PE)*np.cos(psi_PE))
		psi_new = np.arcsin(np.sin(theta_PE)*np.sin(psi_PE)/np.sin(theta_new))+psi
		x = x+steplength*(np.sin(theta_new)*np.cos(psi_new))
		y = y+steplength*(np.sin(theta_new)*np.sin(psi_new))
		z = z+steplength*(np.cos(theta_new))
		theta = theta_new
		psi = psi_new
		return Node(E, x, y, z, theta, psi, steplength)
	else:
		theta_SE = np.arcsin(np.cos(theta_PE))
		E = W
		R_r = np.random.uniform(0,1)
		mfp_el = mfp_elastic(E)
		mfp_in = mfp_inelastic(E)
		mfp_ph = mfp_phonon(E)
		mfp_pol = mfp_polaron(E)
		mfp_T = 1/((1/mfp_el)+(1/mfp_in)+(1/mfp_ph)+(1/mfp_pol))
		steplength = -mfp_T*np.log(R_r)
		psi_SE = np.pi+psi_PE
		theta_new = np.arccos(np.cos(theta)*np.cos(theta_SE)+np.sin(theta)*np.sin(theta_SE)*np.cos(psi_SE))
		psi_new = np.arcsin(np.sin(theta_SE)*np.sin(psi_SE)/np.sin(theta_new))+psi
		x = x+steplength*(np.sin(theta_new)*np.cos(psi_new))
		y = y+steplength*(np.sin(theta_new)*np.sin(psi_new))
		z = z+steplength*(np.cos(theta_new))
		theta = theta_new
		psi = psi_new
		return Node(E, x, y, z, theta, psi, steplength)


def main():
	count_list = []
	E_0 = input1
	theta_0 = input2
	N = input3
	for i in range(N):
		cnt = 0
		E = E_0
		mfp_el = mfp_elastic(E)
		mfp_in = mfp_inelastic(E)
		mfp_ph = mfp_phonon(E)
		mfp_pol = mfp_polaron(E)
		mfp_T = 1/((1/mfp_el)+(1/mfp_in)+(1/mfp_ph)+(1/mfp_pol))
		random.seed(datetime.now())
		R_1 = np.random.uniform(0,1)
		steplength = -mfp_T*np.log(R_1)
		theta = 0
		psi = 0
		theta_PE = theta_0
		psi_PE = 2*np.pi*np.random.uniform(0,1)
		# origin
		x = 0
		y = 0
		z = 0
		theta_new = np.arccos(np.cos(theta)*np.cos(theta_PE)+np.sin(theta)*np.sin(theta_PE)*np.cos(psi_PE))
		psi_new = np.arcsin(np.sin(theta_PE)*np.sin(psi_PE)/np.sin(theta_new))+psi
		x = x+steplength*(np.sin(theta_new)*np.cos(psi_new))
		y = y+steplength*(np.sin(theta_new)*np.sin(psi_new))
		z = z+steplength*(np.cos(theta_new))
		theta = theta_new
		psi = psi_new
		point_1 = Node(E,x,y,z,theta,psi,steplength)
		cnt = buildTree(point_1,100,cnt)
		count_list.append(cnt)
	print 'count_all: ', sum(count_list)
	SEE_yield = sum(count_list)/N
	print 'secondary electron emission yield: ', SEE_yield

if __name__ == "__main__":
	main()
