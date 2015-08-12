#Klassen Model Implementation without Temp. Variance

#Assumptions
#at 300K 
#carrier concentration= Na+Nd


#parameters from tables
######################################
Ymaxp=1410.0
Ymaxb=470.5
Yminp=68.5
Yminb=44.9
Yip=56.1
Yib=29.0
Nref1p=9.2E16
Nref1b=2.23E17
Nref2p=3.41E20
Nref2b=6.1E20
alpha1p=.711
alpha1b=.719
alpha2p=1.98
alpha2b=2.0
s1=.89233
s2=.41372
s3=.19778
s4=.28227
s5=.005978
s6=1.80618
s7=.72169
r1=.7643
r2=2.2999
r3=6.5502
r4=2.3670
r5=-.01552
r6=.6478
########################################

Na=3.50e16 #('Acceptor Concentration =')
Nd=0 #('Donor Concentration =')
Ng=0 #('Gallium Concentration =')
Na=Na+Ng


#Neglected the third term below for Low Dopant Concentrations
#lattice scattering
#Yel=Yminp+((Ymaxp-Yminp)/(1+(pow((Nd/Nref1p),alpha1p))))-(Yip/(1+(pow((Nref2p/Nd),alpha2p))))         
#Yhl=Yminb+((Ymaxb-Yminb)/(1+(pow((Nd/Nref1b),alpha1b))))-(Yib/(1+(pow((Nref2b/Nd),alpha2b))))
Yel=Ymaxp
Yhl=Ymaxb
#Yel=(Yel_1+Yel_2)/2.0
#Yhl=(Yhl_1+Yhl_2)/2.0



kp=(Ymaxp*Ymaxp)/(Ymaxp-Yminp)
kb=(Ymaxb*Ymaxb)/(Ymaxb-Yminb)
mp=(Yminp*Ymaxp)/(Ymaxp-Yminp)
mb=(Yminb*Ymaxb)/(Ymaxb-Yminb)

p=Na-Nd            #mobility hardly changes with changing delta_n
n=0
c=p+n

Ne_sc=Na+Nd+p
Nh_sc=Na+Nd+n


meratio=1.06181427004        # from effective masses program
mhratio=1.32817559494





Pe=(1.36E20)*meratio/c              
Ph=(1.36E20)*mhratio/c

####

def Ge(Pe):
    return (1-(s1/(pow((s2+((pow(meratio,-s4))*Pe)),s3)))+(s5/(pow(Pe,s6))))

def Gh(Ph):
    return (1-(s1/(pow((s2+((pow(mhratio,-s4))*Ph)),s3)))+(s5/(pow(Ph,s6))))

####

me_ratio=meratio/mhratio
mh_ratio=1/me_ratio

####

def Fe(Pe):
    return ((r1*(pow(Pe,r6))+r2+(r3*me_ratio))/((pow(Pe,r6))+r4+(r5*me_ratio)))
def Fh(Ph):
    return ((r1*(pow(Ph,r6))+r2+(r3*mh_ratio))/((pow(Ph,r6))+r4+(r5*mh_ratio)))

####

Ne_sc_eff=Nd+(Ge(Pe)*(Na+(p/Fe(Pe)))
Nh_sc_eff=Na+(Gh(Ph)*Nd)+(n/Fh(Ph))



Ye_D_A_h=float(kp*(Ne_sc/Ne_sc_eff)*((Nref1p/Ne_sc)**alpha1p))+(mp*(c/Ne_sc_eff))
Yh_D_A_e=float(kb*(Nh_sc/Nh_sc_eff)*((Nref1b/Nh_sc)**alpha1b))+(mb*(c/Nh_sc_eff))

#Ye_final=((1/Yel)+(1/Ye_D_A_h)+(1/2000.0))**(-1)               #with correction factor          
Ye_final=((1/Yel)+(1/Ye_D_A_h))**(-1)          

Yh_final=((1/Yhl)+(1/Yh_D_A_e))**(-1)

"""
print "Lattice mobility e: ",Yel
print "Lattice mobility h: ",Yhl

print "Rest e:",Ye_D_A_h
print "Rest h:",Yh_D_A_e
"""
print "\nElectron mobility = ",Ye_final
print "Hole mobility = ",Yh_final

#print 0.73/((Na+Ng-Nd)*1.37*1.6e-19)




