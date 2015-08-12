#Klassen Model Implementation without Temp. Variance

#Assumptions
#at 300K 
#carrier concentration= Na+Nd


#parameters from tables
######################################
Ymaxp=1410
Ymaxb=470.5
Yminp=68.5
Yminb=44.9
Yip=56.1
Yib=29
Nref1p=9.2E16
Nref1b=2.23E17
Nref2p=3.41E20
Nref2b=6.1E20
alpha1p=.711
alpha1b=.719
alpha2p=1.98
alpha2b=2
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

Na=input ('Acceptor Concentration =')
Nd=input ('Donor Concentration =')
Ng=input ('Gallium Concentration =')
#Ng=0
N=Na+Nd+Ng
ND=Na+Ng-Nd

#Neglected the third term below for Low Dopant Concentrations
#lattice scattering
#Yel=Yminp+((Ymaxp-Yminp)/(1+(pow((N/Nref1p),alpha1p))))-(Yip/(1+(pow((Nref2p/N),alpha2p))))         
#Yhl=Yminb+((Ymaxb-Yminb)/(1+(pow((N/Nref1b),alpha1b))))-(Yib/(1+(pow((Nref2b/N),alpha2b))))
Yel=Ymaxp
Yhl=Ymaxb



kp=(Ymaxp*Ymaxp)/(Ymaxp-Yminp)
kb=(Ymaxb*Ymaxb)/(Ymaxb-Yminb)
mp=(Yminp*Ymaxp)/(Ymaxp-Yminp)
mb=(Yminb*Ymaxb)/(Ymaxb-Yminb)

#majority impurity scattering
YeD=(kp*(pow((Nref1p/Nd),alpha1p)))+(mp*ND/Nd)               
YhA=(kb*(pow((Nref1b/Na),alpha1b)))+(mb*ND/Na)


meratio=1.09                    #not sure here (from wikipedia as per my understanding)
mhratio=1.15

Pe=(1.36E20)*meratio/N              #ambiguity here ND or N
Ph=(1.36E20)*mhratio/N

def Ge(Pe):
    return (1-(s1/(pow((s2+((pow(meratio,-s4))*Pe)),s3)))+(s5/(pow(Pe,s6))))

def Gh(Ph):
    return (1-(s1/(pow((s2+((pow(mhratio,-s4))*Ph)),s3)))+(s5/(pow(Ph,s6))))



#minority impurity scattering
YeA=YeD/Ge(Pe)                                       
YhD=YhA/Gh(Ph)

mratio=2   #not sure about this

def Fe(Pe):
    return ((r1*(pow(Pe,r6))+r2+(r3*mratio))/((pow(Pe,r6))+r4+(r5*mratio)))
def Fh(Ph):
    return ((r1*(pow(Ph,r6))+r2+(r3*mratio))/((pow(Ph,r6))+r4+(r5*mratio)))


#electron hole scattering
Yeh=Fe(Pe)*YeD                          
Yhe=Fh(Ph)*YhA


#final values

Yefinal=1/((1/Yel)+(1/YeD)+(1/YeA)+(1/Yeh))
Yhfinal=1/((1/Yhl)+(1/YhD)+(1/YhA)+(1/Yhe))


print ("Electron mobility = ",Yefinal)
print ("Hole mobility = ",Yhfinal)
