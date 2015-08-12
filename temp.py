# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""
Fe Concetration

"""
from math import exp 
import numpy as np


  
z="C:\Users\dell\Desktop\Work\20150610 Ga Fe Flashing\37-148 ref 4 Sam 7 80 Flashes HI.dat"
x="C:\Users\dell\Desktop\Work\20150610 Ga Fe Flashing\37-148 ref 4 Sam 7 before Flashing HI.dat"

###Constants###
Vth=1.1E7
Ev=0.0
Ec=1.16997908699
Kb=1.3806488E-23
T=300
e_ratio_m=1.06181427004
h_ratio_m=1.32817559494
#######


### Input here ###    
Na=input("Boron Conc: " )
#T1=input("Only Fei: ")
#T2=input("Only Feb: ")
#delta_n=input("Excess Conc: ")
def Fe_array():
    delta_n_final=np.genfromtxt(z,usecols=8,skip_header=160,skip_footer=260,delimiter='\t')
    lifetime_final=np.genfromtxt(z,usecols=9,skip_header=160,skip_footer=260,delimiter='\t')
    delta_n_initial=np.genfromtxt(x,usecols=8,skip_header=160,skip_footer=260,delimiter='\t')
    lifetime_initial=np.genfromtxt(x,usecols=9,skip_header=160,skip_footer=260,delimiter='\t')
    
    Fe=[]
    
    for i in range(0,580):
        for j in range(i,580):
            percent_diff= abs(delta_n_final[j]-delta_n_initial[i])/delta_n_initial[i]
            if percent_diff<.1:
                T1=lifetime_final[j]
                T2=lifetime_initial[i]
                delta_n=delta_n_initial[i]
                Fe.append(output(T1,T2,delta_n))
                break
            else:continue
        continue
    return Fe
                
                
######

###Nc and Nv ####
def NC():
    return 2.541E19*(e_ratio_m**1.5)

def NV():
    return 2.541E19*(h_ratio_m**1.5)
    
Nc=NC()
Nv=NV()
#######
    
###Fei and FeB###
class Fei(object):
    
    def __init__(self,delta_n):
        self.sigma_n=1.3E-14
        self.sigma_p=7E-17
        self.delta_e_level=-0.48
        self.delta_n=delta_n
        
    def p1(self):
        return Nv*exp(self.delta_e_level/Kb/T)
    
    def X_Fei(self):
        a=Vth*(Na+self.delta_n)
        b=((Na+self.p1()+self.delta_n)/self.sigma_n)+(self.delta_n/self.sigma_p)
        return a/b
        
class FeB(object):
    
    def __init__(self,delta_n):
        self.sigma_n=5E-15
        self.sigma_p=3E-15
        self.delta_e_level=-0.26
        self.delta_n=delta_n
        
    def n1(self):
        return Nc*exp(self.delta_e_level/Kb/T)
    
    def X_FeB(self):
        a=Vth*(Na+self.delta_n)
        b=((Na+self.delta_n)/self.sigma_n)+((self.n1()+self.delta_n)/self.sigma_p)
        return a/b
    
###########



"""

Space for Auger and other recombination centres

"""

###Output


def output(T1,T2,delta_n):
    g=Fei(delta_n)
    h=FeB(delta_n)
    return (1/(g.X_Fei()-h.X_FeB()))*((1/T1)-(1/T2))

#print "Fei Conc: ",output()
a=Fe_array()
for i in a:
    print i





















