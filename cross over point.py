# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 12:23:52 2015

@author: dell
"""

"""
Cross over point with gallium consideration
over a range of r(ratio of FeB to Fei)
"""

from math import exp
import matplotlib.pyplot as plt


Kb=8.617e-5
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

### objects ###
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

def cop(r,Na):
    a=((r*((FeB.sigma_n*FeB.sigma_p*(FeGa.sigma_n+FeGa.sigma_p))-(FeGa.sigma_n*FeGa.sigma_p*(FeB.sigma_n+FeB.sigma_p))))+(FeGa.sigma_n*FeGa.sigma_p*(FeB.sigma_n+FeB.sigma_p)))    
    
    b=(((Na+FeGa.p1())*r*FeGa.sigma_p*FeB.sigma_n*FeB.sigma_p)+((1-r)*FeGa.sigma_n*FeGa.sigma_p*(FeB.n1()*FeB.sigma_n+Na*FeB.sigma_p)))
         
    c=(Na+FeGa.p1())*(Fei.sigma_p)*((FeB.n1()*FeGa.sigma_n*FeB.sigma_n)+(Na*FeGa.sigma_p*FeB.sigma_p))
    
    cop=(-b+(((b**2.0)+(4.0*a*c))**.5))/2.0/a
    
    return cop
#print a
#print b
#print c
#print "Cross over point:%s"%cop

## plot ##

fig, x1=plt.subplots()
x1.set_xlabel('Ratio of FeB to Fe_Total(dissociated state)')
x1.set_ylabel('Cross over Point (cm-3)')
x1.axis([0,1,1e10,1e15],'normal')
x1.semilogy()

"""plt.switch_backend('wxAgg')"""
for Na,col in zip([1.38e16,1.32e16,9.56e15],['ro','bo','go']):            #using net concentration N_B+N_Ga-N_P
    for r in range (0,11):
        x1.plot(r/10.0,cop(r/10.0,Na),col)
"""
mng=plt.get_current_fig_manager()
mng.frame.Maximize(True)
"""

plt.show()




  
  
     

