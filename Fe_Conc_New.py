# -*- coding: utf-8 -*-
"""
Created on Tue Jul 07 13:59:33 2015

@author: dell
"""

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
import matplotlib.pyplot as plt

#### Files of concern #####
  
after_flashing="C:/Users/dell/Desktop/Work/20150610 Ga Fe Flashing/37-490 ref 4 Sam 7 80 Flashes HI.dat"
before_flashing="C:/Users/dell/Desktop/Work/20150610 Ga Fe Flashing/37-490 ref 4 Sam 7 before Flashing HI.dat"
other="C:/Users/dell/Desktop/Work/20150610 Ga Fe Flashing/37-490 ref 4 Sam 7 60 Flashes HI.dat"

#if initial FeB:Fei ratio is 40:60 then r is 0.4
r=.2

###Constants### ensure floats
Vth=1.1E7
Ev=0.0
Ec=1.16997908699
Kb=8.617e-5     #eV/K
T=300.0
e_ratio_m=1.06181427004
h_ratio_m=1.32817559494
Na=4.9e15    #Net doping

#######


### Input here ###    

# Interpolation
def interpolate(delta_N,L):
    xvals=np.logspace(12,14.3,num=100)      #num represents no. of sampling points in Sample Space
    l_interp=np.interp(xvals,delta_N,L)
    return l_interp

Delta_n=np.logspace(12,14.3,num=100)       # Sample Space

def Fe_array(my_file_1,my_file_2,Fe,l1,l2):
        
    delta_n_final= np.genfromtxt(my_file_1,usecols=8,skip_header=160,skip_footer=550,delimiter='\t')
    lifetime_final= np.genfromtxt(my_file_1,usecols=9,skip_header=160,skip_footer=550,delimiter='\t')
    delta_n_initial= np.genfromtxt(my_file_2,usecols=8,skip_header=160,skip_footer=550,delimiter='\t')
    lifetime_initial= np.genfromtxt(my_file_2,usecols=9,skip_header=160,skip_footer=550,delimiter='\t')
        
    lifetime_final=interpolate(delta_n_final,lifetime_final)    
    lifetime_initial=interpolate(delta_n_initial,lifetime_initial)
    
    for i in range(0,len(lifetime_final)):
        
        T1=lifetime_final[i]
        T2=lifetime_initial[i]
        
        l1.append(T1)                
        l2.append(T2)
        
        Fe.append(abs(output(T1,T2,Delta_n[i])))
        continue
       
    pass
    
    
######

###Nc and Nv for n1 and p1 respectively####
def NC():
    return 2.541E19*(e_ratio_m**1.5)

def NV():
    return 2.541E19*(h_ratio_m**1.5)
    
Nc=NC()
Nv=NV()
#######
    
### Fei and FeB ###
class Fei(object):
    
    def __init__(self,delta_n):
        self.sigma_n=1.3E-14
        self.sigma_p=7E-17
        self.delta_e_level=-0.48
        self.delta_n=delta_n
        
    def p1(self):
        return Nv*exp(self.delta_e_level/Kb/T)
    
    def X_Fei(self):
        a=Vth*(Na+float(self.delta_n))
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
        a=Vth*(Na+float(self.delta_n))
        b=((Na+self.delta_n)/self.sigma_n)+((self.n1()+self.delta_n)/self.sigma_p)
        return a/b
    
###########

###Output###


def output(T1,T2,delta_n):
    g=Fei(delta_n)
    h=FeB(delta_n)
    
    return (1/r/(g.X_Fei()-h.X_FeB()))*((1/T1)-(1/T2))


## cross over point if needed ##

g=Fei(1.1e14)          #independent of delta_n
h=FeB(1.1e14)
cop=(g.sigma_p*Na*((1/h.sigma_n)-(1/g.sigma_n)))+(g.sigma_p/h.sigma_p*h.n1())
#print cop
#print "Fei Conc: ",output()

 
        
#plot individual curves
def plot():
         
    fig, x1=plt.subplots()
    x2=x1.twinx()
    
    Fe_1=[]
    l1_1=[]
    l2_1=[]
    
    Fe_array(after_flashing,before_flashing,Fe_1,l1_1,l2_1)
       
    x1.plot(Delta_n,Fe_1,'go')
    x2.plot(Delta_n,l1_1,'ro')
    x2.plot(Delta_n,l2_1,'bo')
    
    x1.set_xlabel('Delta_n(cm-3) -->')
    x1.set_ylabel('Fe Concentration(cm-3) -->')
    x1.axis([1E11,1E16,1E4,1E14],'normal')
    x1.semilogx()
    x1.semilogy()
    
    x2.set_ylabel('Lifetime(sec) -->')
    x2.semilogy()
    x2.semilogx()
    
    plt.show()  
    pass

#optional for 'other' file       
"""
    Fe_2=[]
    l1_2=[]
    l2_2=[]
    
    Fe_array(other,before_flashing,Fe_2,l1_2,l2_2)
    x1.plot(Delta_n,Fe_2,'bx')
    x2.plot(Delta_n,l1_2,'bo')
    x2.plot(Delta_n,l2_2,'bo')
""" 
    #plt.axvline(cop,color='g')    
    
    
######

plot()

######
