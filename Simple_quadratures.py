import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
from quadratures_Gauss import quadraturaGauss 

def trobaPunts(n,a,b):
    X = np.zeros(n+1)
    for i in range(n+1):
        X[i] = a + (i/(n))*(b-a)
    return X

def trobaPesos(n, a, b, X):
    #Definim la matriu del sistema:
    B = np.ones((n+1,n+1))
    for i in range(1,n+1):
        B[i] = X**i

    #Definim el vector solucio:
    S = np.zeros(n+1)
    for i in range(n+1):
        S[i] = (b**(i+1)-a**(i+1))/(i+1)
    
    #Resolem el sistema per trobar els pesos
    W = np.linalg.solve(B,S)
    return W


def calcula1(w,X):
    integral = 0
    for i in range(X.size):
        integral += w[i]*(np.exp(-(X[i]**2)))
    return integral
  

def calcula2(w,X):
    integral = 0
    for i in range(X.size):
        integral += w[i]*(1/(1+16*X[i]**2))
    return integral  
    
            
    

#Comencem treballant amb les quadratures de Newton-Cotes:
#considerem la integral I1 = sqrt(pi)*erf(1)
#l'aproximarem usant la quadratura N-C amb diferents nombres de punts punts
"""
n = np.arange(1,14)
a = -1
b = 1

Error1 = np.zeros(n.size)
Aprox1 = np.zeros(n.size)
for k in range(n.size): 
    X = trobaPunts(n[k],a,b)
    w = trobaPesos(n[k],a,b,X)
    Aprox1[k] = calcula1(w,X)
    Error1[k] = abs(Aprox1[k] - np.sqrt(np.pi)*scipy.special.erf(1))
    print("Nombre de punts: ", n[k])
    print("L'aproximació numèrica de l'integral és de: ", Aprox1[k])
    print("L'error obtingut és de: ", Error1[k])
    print()
plt.plot(n,np.log10(Error1), '-o', label = 'I1 amb NC')

#Vegem què passa si considerem I2 = 1/2atan(4)

Error2 = np.zeros(n.size)
Aprox2 = np.zeros(n.size)
for k in range(n.size):
    X = trobaPunts(n[k],a,b)
    w = trobaPesos(n[k],a,b,X)
    Aprox2[k] = calcula2(w,X)
    Error2[k] = abs(Aprox2[k] - math.atan(4)/2)
    print("Nombre de punts: ", n[k])
    print("L'aproximació numèrica de l'integral és de: ", Aprox2[k])
    print("L'error obtingut és de: ", Error2[k])
    print()
plt.plot(n,np.log10(Error2), '-o',  label = 'I2 amb NC')

#Tornem a considerar les dues però calcularem les aproximacions usant quadratures de Gauss

Error3 = np.zeros(n.size)
Aprox3 = np.zeros(n.size)
for k in range(n.size):
    X,w = quadraturaGauss(n[k])
    Aprox3[k] = calcula1(w,X)
    Error3[k] = abs(Aprox3[k] - np.sqrt(np.pi)*scipy.special.erf(1))
    print("Nombre de punts: ", n[k])
    print("L'aproximació numèrica de l'integral és de: ", Aprox3[k])
    print("L'error obtingut és de: ", Error3[k])
    print()
plt.plot(n,np.log10(Error3), '-o', label = 'I1 amb Gauss')

Error4 = np.zeros(n.size)
Aprox4 = np.zeros(n.size)
for k in range(n.size):
    X,w = quadraturaGauss(n[k])
    Aprox4[k] = calcula2(w,X)
    Error4[k] = abs(Aprox4[k] - math.atan(4)/2)
    print("Nombre de punts: ", n[k])
    print("L'aproximació numèrica de l'integral és de: ", Aprox4[k])
    print("L'error obtingut és de: ", Error4[k])
    print()
plt.plot(n,np.log10(Error4), '-o', label = 'I2 amb Gauss')
plt.legend(loc='lower left')
"""
