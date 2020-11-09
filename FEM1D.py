# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 11:15:40 2020

author: Luis Carlos A. Rojas Torres
email: luiscarlos.bsf@oceanica.ufrj.br
"""
import numpy as np
import matplotlib.pyplot as plt
from LinearViscoelastic import LinearViscoelastic

def initL(numelem,he):
    L=np.linspace(0,numelem*he,numelem+1)
    return L
def initT(numelem,y):
    T=np.zeros(numelem+1)
    for i in range(T.size):
        T[i]=y
    return T
    
def initK(numelem,he):
    K=np.zeros((numelem,numelem))
    for i in range(numelem):
        K[i][i]=2/he
    for i in range(numelem-1):
        K[i][i+1]=-1/he
        K[i+1][i]=-1/he
    K[0][0]=1/he
    return K    

def fx(PU,lamb,T,f,tau0):
    w=2*np.pi*f
    return 0.5*w*(tau0**2)*PU.get_Jpp(T,f)/lamb

def initF(PU,lamb,T,f,tau0,numelem,he,dy,y):
    F=np.zeros(numelem)
        
    for i in range(numelem-1):
        F[i]+=(he/6)*(2*fx(PU,lamb,T[i],f,tau0)+fx(PU,lamb,T[i+1],f,tau0))
        F[i+1]+=(he/6)*(fx(PU,lamb,T[i],f,tau0)+2*fx(PU,lamb,T[i+1],f,tau0))
        
    F[0]+=dy
    F[-1]+=y/he+((he/6)*(2*fx(PU,lamb,T[-2],f,tau0)+fx(PU,lamb,T[-1],f,tau0)))
    
    return F
    
def submit_static(L,T,K,F,y,numelem,he):
    T[0:-1]=np.linalg.solve(K,F)
    return T
    
