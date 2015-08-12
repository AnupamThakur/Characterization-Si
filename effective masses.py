#Effective masses

Eg_O=1.17
alpha=4.73E-4
beta=636
Eg_T=Eg_O-(((alpha*295)**2)/(295+beta))
e_ratio_m=(6**(2/3.0))*(((0.1905*Eg_O/Eg_T)**2)*0.9163)**(1/3.0)


################
T=300
a=.4435870
b=.3609528E-2
c=.1173515E-3
d=.1263218E-5
e=.3025581E-8
f=.4683382E-2
g=.2286895E-3
h=.7469271E-6
i=.1727481E-8
################

num=(a+(b*T)+(c*T*T)+(d*T*T*T)+(e*T*T*T*T))
den=1+(f*T)+(g*T*T)+(h*T*T*T)+(i*T*T*T*T)

h_ratio_m=((num/den)**2)**(2/3.0)

print Eg_T


print "e_ratio_m: ",e_ratio_m
print "h_ratio_m: ",h_ratio_m



