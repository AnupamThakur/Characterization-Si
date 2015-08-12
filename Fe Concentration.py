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

####Files of concern
  
after_flashing="C:/Users/dell/Desktop/Work/20150610 Ga Fe Flashing/37-232 ref 4 Sam 7 80 Flashes HI.dat"
before_flashing="C:/Users/dell/Desktop/Work/20150610 Ga Fe Flashing/37-232 ref 4 Sam 7 before Flashing HI.dat"
other="C:/Users/dell/Desktop/Work/20150610 Ga Fe Flashing/37-232 ref 4 Sam 7 60 Flashes HI.dat"

###Constants### ensure floats
Vth=1.1E7
Ev=0.0
Ec=1.16997908699
Kb=8.617e-5     #eV/K
T=300.0
e_ratio_m=1.06181427004
h_ratio_m=1.32817559494
Na=6.59e15    #Net doping

#######


### Input here ###    

#T1=input("Only Fei: ")
#T2=input("Only Feb: ")
#delta_n=input("Excess Conc: ")

def Fe_array(my_file_1,my_file_2,Fe,delta_N,l1,l2):
        
    delta_n_final= np.genfromtxt(my_file_1,usecols=8,skip_header=160,skip_footer=260,delimiter='\t')
    lifetime_final= np.genfromtxt(my_file_1,usecols=9,skip_header=160,skip_footer=260,delimiter='\t')
    delta_n_initial= np.genfromtxt(my_file_2,usecols=8,skip_header=160,skip_footer=260,delimiter='\t')
    lifetime_initial= np.genfromtxt(my_file_2,usecols=9,skip_header=160,skip_footer=260,delimiter='\t')
        
        
    for i in range(0,570):
        for j in range(i,570):
            percent_diff= abs(delta_n_final[j]-delta_n_initial[i])/delta_n_initial[i]
            if percent_diff<.01:
                T1=lifetime_final[j]
                T2=lifetime_initial[i]
                delta_n=delta_n_initial[i]
                
                l1.append(T1)                
                l2.append(T2)
                delta_N.append(delta_n)
                Fe.append(abs(output(T1,T2,delta_n)))
                break
            else:continue
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
    
    return (1/(g.X_Fei()-h.X_FeB()))*((1/T1)-(1/T2))


## cross over point if needed ##
"""
g=Fei(1.1e14)
h=FeB(1.1e14)
num=((Na+g.p1())/g.sigma_n)-(Na/h.sigma_n)-(h.n1()/h.sigma_p)
den=(1/h.sigma_n)+(1/h.sigma_p)-(1/g.sigma_n)-(1/g.sigma_p)
print "Cross over point: %s"%(num/den)

#print "Fei Conc: ",output()
"""
# Interpolation 
"""
def interpolate(delta_N,L):
    xvals=np.logspace(1e10,1e16)
    return np.interp(xvals,delta_N,L)
"""        
#plot individual curves
def plot():
         
    fig, x1=plt.subplots()
    x2=x1.twinx()
    
    Fe_1=[]
    delta_N_1=[]
    l1_1=[]
    l2_1=[]
    
    Fe_array(after_flashing,before_flashing,Fe_1,delta_N_1,l1_1,l2_1)
       
    x1.plot(delta_N_1,Fe_1,'go')
    x2.plot(delta_N_1,l1_1,'ro')
    x2.plot(delta_N_1,l2_1,'bo')
    
    x1.set_xlabel('Delta_n(cm-3) -->')
    x1.set_ylabel('Fe_Concentration(cm-3) -->')
    x1.axis([1E11,1E15,1E4,1E14],'normal')
    x1.semilogx()
    x1.semilogy()
    
    x2.set_ylabel('Lifetime(sec) -->')
    x2.semilogy()
    x2.semilogx()
    
#optional for 'other' file       
    """Fe_2=[]
    l1_2=[]
    delta_N_2=[]
    l2_2=[]
    
    Fe_array(other,before_flashing,Fe_2,delta_N_2,l1_2,l2_2)
    x1.plot(delta_N_2,Fe_2,'bo')
    x2.plot(delta_N_2,l1_2,'bo')
    x2.plot(delta_N_2,l2_2,'bo')
"""
######    
    
    plt.show()  
        
    pass

######

plot()

######
