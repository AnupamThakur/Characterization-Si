# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 11:58:51 2015

@author: dell
"""
"""
proves that there is no need for parametrizing compensation
"""

from math import exp
#Finding Po
e_ratio_m=1.06181427004
h_ratio_m=1.32817559494
Eg=1.16997908699
Kb=8.617e-5
T=300.0
Boron=1.07e16
Phospr=4.05e15


def NC():
    return 2.541E19*(e_ratio_m**1.5)

def NV():
    return 2.541E19*(h_ratio_m**1.5)
Nc=NC()
Nv=NV()

Ev=.48          ##Compensation data
Ec=Ev-Eg
def p1():
        return Nv*exp(-Ev/Kb/T)
def n1():
        return Nc*exp(Ec/Kb/T)
n1=n1()
p1=p1()
print n1, p1

Ni=((Nc*Nv)**.5)*exp(-Eg/2/Kb/T)


#alpha is the term exp(-delta_e/Kb/t)
alpha=(-(Nv+4*Phospr)+((((Nv+(4*Phospr))**2)-16*Nv*(Phospr-Boron))**.5))/8/Nv
print alpha
"""
delta_e_ref=log(alpha)*Kb*T

Ea_ratio=1/((Boron/1.3e18)**1.4)
print Ea_ratio

delta_e_new=delta_e_ref*Ea_ratio
"""
