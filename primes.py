#Prime Numbers
a=input("Enter number: ")
b=input("Enter number: ")
primes=[True]*1000000
primes[0]=False
primes[1]=False


i=2
while i*i<b:
    if primes[i]:
        j=i*i
        while j<b:
            primes[j]=False
            j+=i
    i+=1
for i in range(a,b):
    if primes[i]:
        print i
        

