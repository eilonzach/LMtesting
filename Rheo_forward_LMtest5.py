#!/usr/bin/env python
"""This module contains the various subroutines
Usage import """
#######################################################################################

#####################  IMPORT STANDARD MODULES   ######################################   

import sys,os
import numpy as np #for numerical analysis
from scipy import integrate  #integration routine
from subprocess import call #for calling unix commands
from datetime import date  #to give a timestamp to output
import pdb	#for the debugger pdb.set_trace()
from matplotlib.pyplot import * # for the figures
from math import pi,exp,log,sqrt
from ConfigParser import SafeConfigParser
#####################  MODULES   ######################################################   


def PTd_VsQ(freq,P,T,gs,rho,model,highQapprox,Vs0=None,cfile="constants_complete.ini"):
	"""This subroutine calculates the values of Vs and Q at a given 
	PTd condition at a frequency f and density rho. If highQapprox is 1
	approximate values of Vs and Q can be calculated from 
	e.g. McCarthy et al. (2011) eqn. 24, else use their more explicit form eqn. B6"""
	
	global constfileGlob
	constfileGlob = cfile
		
	#### Value of the real and complex part of compliance			
	if Vs0==None:
		J1byJu=getJ1byJu(freq,P,T,gs,model)
		J2byJu=getJ2byJu(freq,P,T,gs,model)
		Ju=getJu(T,P,model)
	else:
		J1byJu=getJ1byJu(freq,P,T,gs,model,Vs0,rho)
		J2byJu=getJ2byJu(freq,P,T,gs,model,Vs0,rho)
		Ju=getJu(T,P,model,Vs0,rho)	

	J1=J1byJu*Ju
	J2=J2byJu*Ju
	
	# EVERYONE AGREES
	qs = J2/J1
	# VS DISAGREEMENT	
	if model=='M11':
		Vs = 1/(sqrt(rho*J1))
	elif model == 'P_M13':
		Vs = 1/(sqrt(rho*J1))
	elif model == 'JF10_eBurg':
		G = pow(J1**2 + J2**2,-0.5)
		Vs = sqrt(G/rho)
	elif model == 'Tak14':
# 		pdb.set_trace()
		G = pow(J1**2 + J2**2,-0.5)
		Vs = sqrt(G/rho)	

	#### When Qs^-1 << 1 (J2/J1 << 1), highQapprox can be 1. If Qs^-1 is not small, Qs^-1 does not equal Qmu^-1 ###### 
	if highQapprox == 0:
		factor=(1+sqrt(1+pow((J2/J1),2)))/2
		qs=qs/factor
		Vs=Vs/sqrt(factor)	

	Q=1/qs
	
	return Vs,Q

def getJ1byJu(freq,P,T,gs,model,Vs0=None,rho=None):
	""" This subroutine gets the value of J1/Ju based on expressions in a paper """
			
#	If model is extended Burgers from Jackson & Faul 2010
	if model == 'JF10_eBurg':
		alpha=getconstant(model,'alpha')
		sigma=getconstant(model,'sigma')
		Delta=getconstant(model,'Delta')	
		tauHR=getconstant(model,'tauHR')	
		tauLR=getconstant(model,'tauLR')	
		tauPR=getconstant(model,'tauPR')
		DeltaP=getconstant(model,'DeltaP')
		ma=getconstant(model,'ma')

		if P > 24.3e9: # LM
			V,Delta,TR,PR,tauF = get_lm_parms()
								
		omega=2*pi*freq	
		tauH=gettau(tauHR,ma,P,T,gs,model)
		tauL=gettau(tauLR,ma,P,T,gs,model)
		tauP=gettau(tauPR,ma,P,T,gs,model)
		
		j1b = (alpha*Delta)/(pow(tauH,alpha) - pow(tauL,alpha))
		j1p = DeltaP/(sigma*sqrt(2*pi))
		i1b = intgrate_j1b(tauL,tauH,alpha,omega,N=1000)
		i1p = integrate.quad(intgrnd_j1p, 0, np.inf, args=(omega,tauP,sigma))[0]

		J1byJu=1.0+(j1b*i1b)+(j1p*i1p)
				
										
	return J1byJu

def getJ2byJu(freq,P,T,gs,model,Vs0=None,rho=None):
	""" This subroutine gets the value of J2 based on expressions in a paper """	
	

#	If model is extended Burgers from Jackson & Faul 2010
	if model == 'JF10_eBurg':
		alpha=getconstant(model,'alpha')
		sigma=getconstant(model,'sigma')
		Delta=getconstant(model,'Delta')	
		tauHR=getconstant(model,'tauHR')	
		tauLR=getconstant(model,'tauLR')	
		tauPR=getconstant(model,'tauPR')
		tauMR=getconstant(model,'tauMR')
		DeltaP=getconstant(model,'DeltaP')
		ma=getconstant(model,'ma')
		mv=getconstant(model,'mv')

		if P > 24.3e9: # LM
			V,Delta,TR,PR,tauF = get_lm_parms()

		omega=2*pi*freq	
		tauH=gettau(tauHR,ma,P,T,gs,model)
		tauL=gettau(tauLR,ma,P,T,gs,model)
		tauP=gettau(tauPR,ma,P,T,gs,model)
		tauM=gettau(tauMR,mv,P,T,gs,model)
		
		j2b = omega*(alpha*Delta)/(pow(tauH,alpha) - pow(tauL,alpha))
		j2p = omega*DeltaP/(sigma*sqrt(2*pi))
		
		tau_arr=np.linspace(tauL,tauH,100)
		i2b = intgrate_j2b(tauL,tauH,alpha,omega,N=1000)
		i2p = integrate.quad(intgrnd_j2p, 0, np.inf, args=(omega,tauP,sigma))[0]
		
		J2byJu = (j2b*i2b)+(j2p*i2p)+(1.0/(omega*tauM))
						
	return J2byJu

def mccarthyXn(taun):
	"""This subroutine gets value of relaxation spectrum at a value of 
	normalized time scale (taun) from McCarthy et al. (2011) eqn. 25 """
	if taun >= 1.0e-11:
		Xn=0.32*(pow(taun,(0.39-0.28/(1+2.6*pow(taun,0.1)))))
	else:
		Xn=1853.0*sqrt(taun)
	return Xn
			
def getfn(freq,P,T,gs,model,Vs0 = None,rho=None):
	"""Get the normalized frequency from McCarthy et al. (2011) eqn. 19 and 22"""
	eta0=getconstant(model,'eta0')
	m=getconstant(model,'m')
	
	if Vs0==None:
		Ju=getJu(T,P,model)
	else:
		Ju=getJu(T,P,model,Vs0,rho)	
		
	tau0=Ju*eta0
	taur=gettau(tau0,m,P,T,gs,model)
	if model == "P_M13":
		taur = gettauM_P_M13(T,P,Ju)
	
	fn=freq*taur	
	return fn

	
def getJu(T,P,model,Vs0 = None,rho = None ):
	"""This gets the unrelaxed modulus at T,P conditions based on constants from study"""

	if Vs0 is None:
		GUR=getconstant(model,'GUR')
		dGdT=getconstant(model,'dGdT') 
		dGdP=getconstant(model,'dGdP')
		TR=getconstant(model,'TR')
		PR=getconstant(model,'PR')
		Ju=1/(GUR + dGdT*(T-TR) + dGdP*(P-PR))
	else:
		Ju=1/(rho*pow(Vs0,2))	
	
	return Ju

def gettau(tau0,m,P,T,gs,model):
	"""calculate the relaxation time at P,T,gs, given model and reference(0) relaxation 
	time and grainsize exponent"""	
	PR=getconstant(model,'PR')
	TR=getconstant(model,'TR') 
	gsR=getconstant(model,'gsR')
	E=getconstant(model,'E')
	V=getconstant(model,'V')
	R=getconstant('Common','R')	
	
	if P > 24.3e9: # LM
		V,Delta,TR,PR,tauF = get_lm_parms()
		
# 	print V
# 	print Delta
# 	print TR
# 	print PR
# 	print tauMR
# 	pdb.set_trace()
		
	tau=tau0*tauF*pow((gs/gsR),m)*exp( (E/R)*(TR-T)/(TR*T) + (V/R)*(P*TR-PR*T)/(T*TR) )
	
	return tau

		
#####  Integrations for JF10  #####
def intgrate_j1b(limL,limH,alpha,omega,N=1000):
	xx = np.logspace(np.log10(limL),np.log10(limH),N)
	yy = np.zeros(len(xx))
	for ii in range(0,len(xx)):
		yy[ii] = pow(xx[ii],alpha-1)/(1 + pow(omega*xx[ii],2))
		
	I=integrate.cumtrapz(yy,xx)
	I = I[len(I)-1]
	return I
	
def intgrate_j2b(limL,limH,alpha,omega,N=1000):
	xx = np.logspace(np.log10(limL),np.log10(limH),N)
	yy = np.zeros(len(xx))
	for ii in range(0,len(xx)):
		yy[ii] = pow(xx[ii],alpha)/(1 + pow(omega*xx[ii],2))
		
	I=integrate.cumtrapz(yy,xx)
	I = I[len(I)-1]
	return I
	
def intgrate_j1p(limL,limH,omega,tauP,sigma,N=1000):
	xx = np.logspace(np.log10(limL),np.log10(limH),N)
	yy = np.zeros(len(xx))
	for ii in range(0,len(xx)):
		yy[ii] = (1/xx[ii])*(1/(1 + pow(omega*xx[ii],2)))*exp(-pow(log(xx[ii]/tauP),2)/(2*pow(sigma,2)))
		
	I=integrate.cumtrapz(yy,xx)
	I = I[len(I)-1]
	return I
	
		
def intgrnd_j1b(tau,alpha,omega):
	return pow(tau,alpha-1)/(1. + pow(omega*tau,2.))

def intgrnd_j1p(tau,omega,tauP,sigma):
	return (1/tau)*(1/(1 + pow(omega*tau,2)))*exp(-pow(log(tau/tauP),2)/(2*pow(sigma,2)))

def intgrnd_j2b(tau,alpha,omega):
	return pow(tau,alpha)/(1. + pow(omega*tau,2.))
	
def intgrnd_j2p(tau,omega,tauP,sigma):
	return (1/(1 + pow(omega*tau,2)))*exp(-pow(log(tau/tauP),2)/(2*pow(sigma,2)))



def getconstant(model,par):
	"""This subroutine gets the value of a constant relevant for a modeling 
	approach"""
	parser = SafeConfigParser()
# 	parser.read('constants_complete_LMtest4.ini')
	parser.read(constfileGlob)
	val=float(parser.get(model,par))
	return val


def get_lm_parms():
	V=getconstant('LOWM','lmV')
	Delta=getconstant('LOWM','lmDelta')
	TR=getconstant('LOWM','lmTR')
	PR=getconstant('LOWM','lmPR')
	tauF=getconstant('LOWM','lmtauF')
	return V,Delta,TR,PR,tauF


	
#########################  	PLOTTING   ############################################   	
