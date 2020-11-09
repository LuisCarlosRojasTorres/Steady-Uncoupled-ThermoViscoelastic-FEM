# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 13:06:22 2020

author: Luis Carlos A. Rojas Torres
email: luiscarlos.bsf@oceanica.ufrj.br
"""

import numpy as np

class LinearViscoelastic:
    '''DOCUMENTATION: A library for a linear viscoelastic solid represented 
    by a n-th order Prony Serie'''

    def __init__(self,stress,rTimes,v=0.45,Tref=23,c1=7.32,c2=93.62):
        self.stress = np.array(stress)*10**6   #Array of Elastic Modules
        self.rTimesTref = np.array(rTimes) ##Array of Relaxation Times at Tref
        self.rTimes = np.array(rTimes)     #Array of Relaxation Times at T
        self.T = Tref
        self.Tref = Tref              #Reference temperature
        self.c1 = c1                  #WLF coefficient c1
        self.c2 = c2                  #WLF coefficient c2
        self.v = v                    #Poisson's Ratio 

    def modify_T(self,T):
        #Change the temperature of analysis
	    self.T=T                    #Given temperature
	    self.rTimes=self.rTimesTref*self.get_shiftFactor()
	
    def get_shiftFactor(self):
        #Calculates the shift factor WLF at a given temperature (T)
    	return 10**(-self.c1*(self.T-self.Tref)/(self.c2+self.T-self.Tref))

    def get_Er(self,T,time=0):
        Er=self.stress[0]
        
        self.modify_T(T)
        for i in range(1,len(self.stress)):
            Er=Er+self.stress[i]*np.exp(-time/self.rTimes[i-1])
        return Er    

    def get_Ep(self,T,frec=1):
        w=2*np.pi*frec
        Einf=self.stress[0]
        Ep=0
        self.modify_T(T)
        for i in range(1,len(self.stress)):
            Ep+=self.stress[i]*(w*self.rTimes[i-1])**2/(1+(w*self.rTimes[i-1])**2)
        return Einf+Ep

    def get_Epp(self,T,frec=1):
        w=2*np.pi*frec
        Epp=0
        self.modify_T(T)
        for i in range(1,len(self.stress)):
            Epp+=(self.stress[i]*w*self.rTimes[i-1])/(1+(w*self.rTimes[i-1])**2)
        return Epp
    
    def get_Dp(self,T,frec):
        Dp=self.get_Ep(T,frec)
        Dp/=(self.get_Ep(T,frec)**2+self.get_Epp(T,frec)**2)
        return Dp
    
    def get_Dpp(self,T,frec):
        Dpp=self.get_Epp(T,frec)
        Dpp/=(self.get_Ep(T,frec)**2+self.get_Epp(T,frec)**2)
        return Dpp
    
    def get_Jpp(self,T,frec):
        Jpp=2*((1+self.v)*self.get_Epp(T,frec))
        Jpp/=(self.get_Ep(T,frec)**2+self.get_Epp(T,frec)**2)
        return Jpp
        