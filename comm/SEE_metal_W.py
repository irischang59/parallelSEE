# -*- coding: utf-8 -*-
# Discrete Monte Carlo Model

from __future__ import division
import sys
import numpy as np
import random
from datetime import datetime
# from sympy.solvers import solve
# from sympy import Symbol


if len(sys.argv) != 4:
	print "The Single Scattering Monte Carlo Model"
	exit("Usage: %s [primary electron energy (keV)] [incident angle (radian)] [sampling times] " % sys.argv[0])

input1 = float(sys.argv[1])              # primary electron energy [eV]
input2 = float(sys.argv[2])              # incident angle [radian]
input3 = int(sys.argv[3])                # sampling times


'''
########## Copper (Cu) ##########
Z = 29.0                                 # atomic number
Na = 6.022*10**23                        # Avagadro's number
rho = 8.92                               # density [g/cm^3]
AW = 63.546                              # atomic weight [g/mol]
N = 8.5*10**22                           # volume number density [number of atoms/cm^3]
Nv = 11
E_F = 7*10**(-3)                         # Fermi energy [keV]
phi = 4.45*10**(-3)                      # work function [keV]
V = E_F+phi                              # surface potential barrier
E_B =                                    # binding energy [keV]
a1 = -2.11246
b1 = -1.1456
c1 = 0.22332
d1 = -0.620971
e = 2.7813
a2 = 0.712004
b2 = 0.00196709
c2 = -0.0218025
d2 = -0.059572
'''

'''
########## Gold (Au) ##########
Z = 79.0                                 # atomic number
Na = 6.022*10**23                        # Avagadro's number
rho = 19.3                               # density [g/cm^3]
AW = 196.97                              # atomic weight [g/mol]
N = 6.3*10**23                           # volume number density [number of atoms/cm^3]
Nv = 11
E_F = 5.5*10**(-3)                       # Fermi energy [keV]
phi = 4.7*10**(-3)                       # work function [keV]
V = E_F+phi                              # surface potential barrier
E_B =                                    # binding energy [keV]
a1 = -2.55333
b1 = -0.701946
c1 = 0.0717782
d1 = -0.197945
e = 2.7813
a2 = 0.711771
b2 = 0.00193186
c2 = -0.0260514
d2 = -0.0609567
'''

########## Tungsten (W) ##########
Z = 74.0                                 # atomic number
Na = 6.022*10**23                        # Avagadro's number
rho = 19.25                              # density [g/cm^3]
AW = 183.84                              # atomic weight [g/mol]
N = 6.31*10**22                          # volume number density [number of atoms/cm^3]
Nv = 6.0                                 # number of valence electrons per atom or molecule
z = 2.0                                  # electron number density [number of electrons/cm^3]
Ns = 6.0
n = N*Ns
k_F = 6.66511*10**7*(n/(10**22))*(1.0/3.0)   # wave number [cm^-1]
E_F = 3.64645*10**(-15)*(n**(2.0/3.0))   # Fermi energy [eV]
E_F = E_F*10**(-3)                       # Fermi energy [keV]
# print 'E_F: ', E_F
phi = 4.55*10**(-3)                      # work function [keV]
V = E_F+phi                              # surface potential barrier
# print 'V: ', V*10**3
E_B = 10.2                               # binding energy [keV]
a1 = -2.0205
b1 = -1.2589
c1 = 0.271737
d1 = -0.695477
e = 2.7813
a2 = 0.71392
b2 = 0.00197916
c2 = -0.0172852
d2 = -0.0570799


########## Physical Constants ##########
r_Bohr = 5.29*10**(-9)                   # Bohr radius [cm]
e = 1.6*10**(-19)                        # electron charge [C]
epsilon = 8.85*10**(-12)                 # permittivity of vaccum [F/m]
m_e = 9.1*10**(-31)                      # mass of electron [kg]
omega_p = 56.4*np.sqrt(n*10**6)          # plasma frequency = np.sqrt(n*e**2/(epsilon*m_e)) [rad/s]
h_bar = 6.58*10**(-16)                   # reduced Planck constant [eV-s/rad]
E_plasmon = h_bar*omega_p                # plasmon energy [eV]


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
	# print data
	if height == 0:
		# print "Terminated."
		return count
	elif root.E <= V:
		# print "The energy of secondry electron falls below surface potential barrier."
		return count
	elif root.z < 0:
		# print "Secondry electron emitted out of the metal surface."
		count_out = count+1
		count = count_out
		data = np.asarray([1000*(root.E-V), root.x, root.z, root.theta]).reshape(1,4)
		print data
		return count

	mfp_el = mfp_elastic(root.E)
	# print "elsetic mean free path: ", mfp_el
	mfp_p = mfp_plasmon(root.E)
	# print 'mfp_p', mfp_p
	mfp_c = mfp_conduction(root.E)
	# print 'mfp_c', mfp_c
	# mfp_core = mfp_Gr_M(root.E)
	# mfp_core = mfp_Gr_K(root.E)+mfp_Gr_L(root.E)+mfp_Gr_M(root.E)
	# print 'mfp_shell', mfp_shell
	# mfp_T = 1/((1/mfp_el)+(1/mfp_p)+(1/mfp_c)+(1/mfp_core))
	mfp_T = 1/((1/mfp_el)+(1/mfp_p)+(1/mfp_c))
	# print "total mean free path: ", mfp_T
	Test_scatter = np.random.uniform(0,1) 
	P_el = mfp_T/mfp_el
	# print 'P_el', P_el
	if Test_scatter > P_el:
		# print 'inelastic scattering'
		R_E = np.random.uniform(0,1)                    # random number
		if Test_scatter < (1/mfp_el+1/mfp_p)/(1/mfp_T):
			# print 'plasmon excitation'
			next = plasmonExcitation(root, root.E, root.x, root.y, root.z, root.theta, root.psi, root.steplength)
			root.left = next
			root.right = None
			count_left = buildTree(next, height-1, count)
			count = count_left
			return count
		# elif Test_scatter < (1/mfp_el+1/mfp_p+1/mfp_c)/(1/mfp_T):
		else:
			# print 'conduction electron excitation'
			left = generateSE_conduction(root, root.E, root.x, root.y, root.z, root.theta, root.psi, root.steplength, True, R_E)
			right = generateSE_conduction(root, root.E, root.x, root.y, root.z, root.theta, root.psi, root.steplength, False, R_E)
			root.left = left
			root.right = right
			count_left = buildTree(left, height-1, count)
			count = count_left
			count_right = buildTree(right, height-1, count)
			count_leftandright = count_right
			count = count_leftandright
			return count
			'''
		else:
			# print 'core-shell electron excitation'
			left = generateSE_core(root, root.E, root.x, root.y, root.z, root.cx, root.cy, root.cz, root.theta, root.psi, root.steplength, True, R_E)
			right = generateSE_core(root, root.E, root.x, root.y, root.z, root.cx, root.cy, root.cz, root.theta, root.psi, root.steplength, False, R_E)
			root.left = left
			root.right = right
			count_left = buildTree(left, height-1, count)
			count = count_left
			count_right = buildTree(right, height-1, count)
			count_leftandright = count_right
			count = count_leftandright
			return count
			'''
	else:
		# print 'elastic scattering'
		next = elasticScatter(root, root.E, root.x, root.y, root.z, root.theta, root.psi, root.steplength)
		root.left = next
		root.right = None
		count_left = buildTree(next, height-1, count)
		count = count_left
		return count

def mfp_elastic(E):
	sigma_el = 3*10**(-18)*(Z**1.7)/(E+0.005*(Z**1.7)*np.sqrt(E)+0.0007*Z**2/np.sqrt(E)) 
	return AW/(Na*rho*sigma_el)*10**8   

def mfp_conduction(E):
	Na = 6.022*10**23                            # Avagadro's number
	rho = 19.25                                  # bulk density [g/cm^3]
	AW = 183.84                                  # atomic weight [g/mol]
	N = 6.3*10**22                               # volume number density [number of atoms/cm^3]
	z = 2.0
	Ns = 6.0
	n = N*Ns
	E_F = 3.64645*10**(-15)*(n**(2.0/3.0))       # Fermi energy [eV]
	phi = 4.55                                   # work function [eV]
	sigma_c = AW/(rho*Na)*6.5141*10**(-14)*Na/(1000*E)*((1000*E)-E_F-phi)/(phi*((1000*E)-E_F))     # conduction electron cross section [cm^2]
	return AW/(Na*rho*sigma_c)*10**8             # mean free path for conduction electron [A]

def mfp_plasmon(E):
	theta_p = E_plasmon/(2*1000*E)
	sigma_p = z*AW*theta_p/(2*Na*rho*r_Bohr)*(np.log(theta_p**2+0.175**2)-np.log(theta_p**2))
	return AW/(Na*rho*sigma_p)*10**8 
	# return (2*(r_Bohr*10**10)*(1000*E)/E_plasmon)*np.log((np.sqrt(1000*E_F+E_plasmon)-np.sqrt(1000*E_F))/(np.sqrt(1000*E)-np.sqrt(1000*E-E_plasmon)))

def mfp_Gr_K(E):
	N_K = 2.0
	Eb_K = 69525.0
	t = 1000*E
	U = t/Eb_K
	sigma_Gr_K = 6.5141*10**(-14)*N_K*(1/Eb_K**2)*(1/U)*((U-1)/(U+1))**(3.0/2.0)*(1+2.0/3.0*(1-1/(2*U))*np.log(2.7+np.sqrt(U-1)))
	return AW/(Na*rho*sigma_Gr_K)*10**8 

def mfp_Gr_L(E):
	N_L = 8.0
	Eb_L = 12100.0
	t = 1000*E
	U = t/Eb_L
	sigma_Gr_L = 6.5141*10**(-14)*N_L*(1/Eb_L**2)*(1/U)*((U-1)/(U+1))**(3.0/2.0)*(1+2.0/3.0*(1-1/(2*U))*np.log(2.7+np.sqrt(U-1)))
	return AW/(Na*rho*sigma_Gr_L)*10**8

def mfp_Gr_M(E):
	N_M = 18.0
	Eb_M = 2820.0
	t = 1000*E
	U = t/Eb_M
	sigma_Gr_M = 6.5141*10**(-14)*N_M*(1/Eb_M**2)*(1/U)*((U-1)/(U+1))**(3.0/2.0)*(1+2.0/3.0*(1-1/(2*U))*np.log(2.7+np.sqrt(U-1)))
	return AW/(Na*rho*sigma_Gr_M)*10**8

def mfp_shell(E):
	return mfp_Gr_K(E)+mfp_Gr_L(E)+mfp_Gr_M(E)

def elasticScatter(root, E, x, y, z, theta, psi, steplength): 
	mfp_el = mfp_elastic(E)                               # elastic mean free path [A]
	mfp_p = mfp_plasmon(E)
	mfp_c = mfp_conduction(E)
	mfp_T = 1/((1/mfp_el)+(1/mfp_p)+(1/mfp_c))
	R1 = np.random.uniform(0,1)
	steplength = -mfp_T*np.log(R1)
	alpha = 10**(a1+b1*np.log10(E)+c1*np.log10(E)**2+d1/(e**np.log10(E)))   # screening factor
	beta_star = a2+b2*np.sqrt(E)*np.log(E)+c2*np.log(E)/E+d2/E
	if beta_star > 1:
		beta = 1
	else:
		beta = beta_star
	R_2 = np.random.uniform(0,1)                         # random number
	R_max = (np.cos(np.deg2rad(180**beta))+alpha*np.cos(np.deg2rad(180**beta))-1-alpha)/(np.cos(np.deg2rad(180**beta))-1-2*alpha)
	R_star = R_2*R_max
	theta_PE = np.deg2rad(np.rad2deg(np.arccos(1-(2*alpha*R_star)/(1+alpha-R_star)))**(1/beta))
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


def plasmonExcitation(root, E, x, y, z, theta, psi, steplength):
	E = E-(E_plasmon/1000)
	mfp_el = mfp_elastic(E)                               # elastic mean free path [A]
	mfp_p = mfp_plasmon(E)
	mfp_c = mfp_conduction(E)
	mfp_T = 1/((1/mfp_el)+(1/mfp_p)+(1/mfp_c))
	R_r = np.random.uniform(0,1)           # random number
	steplength = -mfp_T*np.log(R_r)
	x = x+steplength*np.sin(theta)*np.cos(psi)
	y = y+steplength*np.sin(theta)*np.sin(psi)
	z = z+steplength*np.cos(theta)
	return Node(E, x, y, z, theta, psi, steplength)


def generateSE_conduction(root, E, x, y, z, theta, psi, steplength, isGenerateLeft, R_E): 
	R_s = np.random.uniform(0,1)                       # random number
	A = (E-E_F)/(E-E_F-phi)
	R_3 = np.random.uniform(0,1)                       # random number
	dE = (R_s*E_F-A*(E_F+phi))/(R_s-A)
	R = random.uniform(-1,1)
	if R>=0:
		theta_PE = np.arcsin(np.sqrt(dE/E))
	else:
		theta_PE = np.pi-np.arcsin(np.sqrt(dE/E))
	psi_PE = 2*np.pi*R_3
	R = np.random.uniform(0,1) 
	if isGenerateLeft:
		E = E-dE
		mfp_el = mfp_elastic(E)                               # elastic mean free path [A]
		mfp_p = mfp_plasmon(E)
		mfp_c = mfp_conduction(E)
		mfp_T = 1/((1/mfp_el)+(1/mfp_p)+(1/mfp_c))
		R_r = np.random.uniform(0,1)                   # random number
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
		E = dE
		mfp_el = mfp_elastic(E)                               # elastic mean free path [A]
		mfp_p = mfp_plasmon(E)
		mfp_c = mfp_conduction(E)
		mfp_T = 1/((1/mfp_el)+(1/mfp_p)+(1/mfp_c))
		R_r = np.random.uniform(0,1)           # random number
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

'''
def generateSE_core(root, E, x, y, z, theta, psi, steplength, isGenerateLeft, R_E): 
	R_c = np.random.uniform(0,1)
	R_3 = np.random.uniform(0,1)                   # random number
	dE = Symbol('dE')
	solve((E*(1-2*E_B/(dE+E_B)**(3/2))*(1+(2/3-1/3*E_B/dE)*np.log(2.7+np.sqrt(dE/E_B-1))))/(dE*(1-2*E_B/(E+E_B)**(3/2))*(1+(2/3-1/3*E_B/E)*np.log(2.7+np.sqrt(dE/E_B-1))))-R_c,dE)
	theta_PE = np.arcsin(np.sqrt(dE/E))
	if isGenerateLeft:
		theta = theta_PE
		E = E-dE
		mfp = mfp_Gr_M(E)
		R_r = np.random.uniform(0,1)                   # random number
		steplength = -mfp*np.log(R_r)
		psi = 2*np.pi*R_3
		AM = -cx/cz
		AN = 1/(np.sqrt(1+AM**2))
		V1 = AN*np.sin(theta)
		V2 = AM*AN*np.sin(theta)
		V3 = np.cos(psi)
		V4 = np.sin(psi)
		ca = (cx*np.cos(theta))+(V1*V3)+(cy*V2*V4)
		cb = (cy*np.cos(theta))+(V4*(cz*V1-cx*V2))
		cc = (cz*np.cos(theta))+(V2*V3)-(cy*V1*V4)
		x = x+steplength*ca
		y = y+steplength*cb
		z = z+steplength*cc
		cx = ca
		cy = cb
		cz = cc
		theta = np.arcsin(cx)
		return Node(E, x, y, z, theta, psi, steplength)
	else:
		theta = np.arcsin(np.cos(theta_PE))
		E = dE
		mfp = mfp_Gr_M(E)
		R_r = np.random.uniform(0,1)                   # random number
		steplength = -mfp*np.log(R_r)
		psi = np.pi+2*np.pi*R_3
		AM = -cx/cz
		AN = 1/(np.sqrt(1+AM**2))
		V1 = AN*np.sin(theta)
		V2 = AM*AN*np.sin(theta)
		V3 = np.cos(psi)
		V4 = np.sin(psi)
		ca = (cx*np.cos(theta))+(V1*V3)+(cy*V2*V4)
		cb = (cy*np.cos(theta))+(V4*(cz*V1-cx*V2))
		cc = (cz*np.cos(theta))+(V2*V3)-(cy*V1*V4)
		x = x+steplength*ca
		y = y+steplength*cb
		z = z+steplength*cc
		cx = ca
		cy = cb
		cz = cc
		theta = np.arcsin(cx)
		return Node(E, x, y, z, theta, psi, steplength)
'''

def main():
	count_list = []
	E_0 = input1 #("The energy of the primary electron beam: ")
	theta_0 = input2 #("The incident angle of the primary electron beam: ")
	N = input3 #("Sampling times: ")
	for i in range(N):
		cnt = 0
		E = E_0
		mfp_el = mfp_elastic(E)
		mfp_p = mfp_plasmon(E)
		mfp_c = mfp_conduction(E)
		mfp_T = 1/((1/mfp_el)+(1/mfp_p)+(1/mfp_c))
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



