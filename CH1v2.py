# -*- coding: utf-8 -*-
"""
STEADY STATE
Resolve através do Metodo dos Elementos Finitos o problema
termo-vicoelastico de geracao de temperatura numa laje submetida a 
carregamento cizalhante.
author: Luis Carlos A. Rojas Torres
email: luiscarlos.bsf@oceanica.ufrj.br
"""
import numpy as np
import FEM1D as f1
from LinearViscoelastic import LinearViscoelastic

XFrec=list()
YTemp=list()
for i in range(1,19):
    f=1+0.5*i
    tau0=100000
    lamb=0.214
    ite=10
    print("Frequency: ",f,"Hz")
    elastic_modulus = [38,20.14,14.3,9.97,7.25,5.15,3.4,2.1,1.8] #in MPa
    relaxation_times =[0.58,3.13,18.72,125.41,1042,10942,159569,5215397]    #in seconds

    PU=LinearViscoelastic(elastic_modulus,relaxation_times)
    lenght=0.1
    numelem = 20
    he=lenght/numelem
    dy=0
    y=23

    L=f1.initL(numelem,he)
    T=f1.initT(numelem,y)
    K=f1.initK(numelem,he)

    for num_ite in range(1,ite+1):
        F=f1.initF(PU,lamb,T,f,tau0,numelem,he,dy,y)
        T=f1.submit_static(L,T,K,F,y,numelem,he)
        
    XFrec.append(f)
    YTemp.append(T[0])

