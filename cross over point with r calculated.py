# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 10:40:17 2015

@author: dell
"""
# -*- coding: utf-8 -*-
"""
Cross over point with gallium consideration

not for samples without gallium

r calculated from boltzmann statistics assuming thermal equilibrium  (1.4830227 by Bartel paper)

This is explicitly for the case when sample Fe-A is dissociated first and then allowed to rest (acc. to me)

"""
from math import exp


Kb=8.617e-5         # eV/K
T=300.0
Ec=1.16997908699
e_ratio_m=1.06181427004
h_ratio_m=1.32817559494
###Nc and Nv for n1 and p1 respectively####
def NC():
    return 2.541E19*(e_ratio_m**1.5)

def NV():
    return 2.541E19*(h_ratio_m**1.5)
    
Nc=NC()
Nv=NV()
#######

### 'Defects' objects ###
class Fei(object):
    
    def __init__(self):
        self.sigma_n=1.3E-14
        self.sigma_p=7E-17
        self.delta_e_level=-0.48
        
                    
class FeB(object):
    
    def __init__(self):
        self.sigma_n=5E-15
        self.sigma_p=3E-15
        self.delta_e_level=-0.26
        
    def n1(self):
        return Nc*exp(self.delta_e_level/Kb/T)
    
class FeGa(object):
    
    def __init__(self):
        self.sigma_n=4e-14
        self.sigma_p=2e-14
        self.delta_e_level=-.2
    
    def p1(self):
        return Nv*exp(self.delta_e_level/Kb/T)
    

Fei=Fei()
FeB=FeB()
FeGa=FeGa()
    
## For quadratic polynomial : ax2+bx-c=0 ##
# Here

#Ratio of FeB to Fei is r
N_B=3.23e16
N_Ga=8.81e15
delta_e=0.029   #can be changed.... this is binding energy difference at 300 K according to paper
Na=1.32e16              #use net conc. here  N_B+N_Ga-N_P

def r():
    z=(N_B*exp(-delta_e/Kb/T))+N_Ga
    return N_B*exp(-delta_e/Kb/T)/z

r=r()    

def cop():
    a=((r*((FeB.sigma_n*FeB.sigma_p*(FeGa.sigma_n+FeGa.sigma_p))-(FeGa.sigma_n*FeGa.sigma_p*(FeB.sigma_n+FeB.sigma_p))))+(FeGa.sigma_n*FeGa.sigma_p*(FeB.sigma_n+FeB.sigma_p)))    
    
    b=(((Na+FeGa.p1())*r*FeGa.sigma_p*FeB.sigma_n*FeB.sigma_p)+((1-r)*FeGa.sigma_n*FeGa.sigma_p*(FeB.n1()*FeB.sigma_n+Na*FeB.sigma_p)))
         
    c=(Na+FeGa.p1())*(Fei.sigma_p)*((FeB.n1()*FeGa.sigma_n*FeB.sigma_n)+(Na*FeGa.sigma_p*FeB.sigma_p))
    
    cop=(-b+(((b**2.0)+(4.0*a*c))**.5))/2.0/a
    
    return cop
#print a
#print b
#print c
print r
print "Cross over point: %s"%cop()
