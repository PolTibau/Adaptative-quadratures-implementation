import numpy as np
from quadratures_simples import trobaPesos
from quadratures_simples import trobaPunts
from quadratures_compostes import SimpsonCompost
from quadratures_compostes import trapeziCompost
from quadratures_compostes import calcula

import matplotlib.pyplot as plt
import scipy

def f(x):
    return np.sin(np.exp(2*x))

#EXERCICI 1
#Volem calcular la integral de f(x) entre a = 0 i b = 2 amb 6 xifres signif.
a = 0
b = 2

print("La integral resultant es correspon a: ", scipy.integrate.quad(f,a,b)[0])
#Comencem calculant una aproximacio amb Simpson i trapezi compostos
M = np.array([4,8,16,32,64,128])
T = np.zeros(M.size)
S = np.zeros(M.size)
ET = []
ES = []
for i in range(M.size):
    T[i] = trapeziCompost(f, a, b, M[i])
    ET.append(abs(T[i] - scipy.integrate.quad(f,a,b)[0]))
    S[i] = SimpsonCompost(f, a, b, M[i])
    ES.append(abs(S[i] - scipy.integrate.quad(f,a,b)[0]))
    print("L'aproximació pel mètode del trapezi amb ", M[i], "subintervals és ", T[i])
    print("L'aproximació pel mètode de Simpson amb ", M[i], "suvintervals és ", S[i])
    print()
print(ET[-1], ES[-1])
plt.title('Evolucio errors amb Trapezi i Simpson')
plt.plot(np.log(M),np.log10(ET), '-o', label = 'Ev error Trapezi')
plt.plot(np.log(M),np.log10(ES), '-o', label = 'Ev error Simpson')
plt.legend()
plt.show()


#EXERCICI 2
AproxT = False
AproxS = False
compteT = 32
compteS = 32

while(not AproxT):
    t = trapeziCompost(f, a, b, compteT)
    if(abs(t - scipy.integrate.quad(f,a,b)[0]) < 1e-6):
        AproxT = True
    else:
        compteT += 1
    
while(not AproxS):
    s = SimpsonCompost(f, a, b, compteS)
    if(abs(s - scipy.integrate.quad(f,a,b)[0]) < 1e-3):
        AproxS = True
    else:
        compteS += 1 
    
print("Es necessiten", compteT, "subintervals per aconseguir 6 xifres significatives amb el mètode del Trapezi")
print("Es necessiten", compteS, "subintervals per aconseguir 6 xifres significatives amb el mètode de Simpson")
print("L'aproximació pel mètode del trapezi amb ", compteT, "subintervals és ", t)
print("L'aproximació pel mètode de Simpson amb ", compteS, "subintervals és ", s)
print()

#EXERCICI 4
#Plantegem una quadratura adaptativa,
def SimpsonAdaptativa(f,a,b,eps):
    I = SimpsonCompost(f, a, b, 1)
    II = SimpsonCompost(f, a, (a+b)/2, 1)
    III = SimpsonCompost(f, (a+b)/2, b, 1)
    Eab = abs(I - (II + III))
    if(Eab < eps*(b-a)):
        return I
    else:
        return SimpsonAdaptativa(f, a, (a+b)/2,eps) + SimpsonAdaptativa(f, (a+b)/2, b, eps)


tol1 = 0.5e-3
tol2 = 0.5e-6
AP1 = SimpsonAdaptativa(f, a, b, tol1)
AP2 = SimpsonAdaptativa(f, a, b, tol2)
print("Usant un mètode recursiu i una tolerància de 1e-3 obtenim ", AP1)
print("Usant un mètode recursiu i una tolerància de 1e-6 obtenim ", AP2)
print("El primer resultat té un error absolut de: " , abs(AP1 - scipy.integrate.quad(f,a,b)[0]))
print("El segon resultat té un error absolut de: " , abs(AP2 - scipy.integrate.quad(f,a,b)[0]))

#EXERCICI 5
def SimpsonAdaptativa2(f,a,b,eps,ABS):
    I = SimpsonCompost(f, a, b, 1)
    II = SimpsonCompost(f, a, (a+b)/2, 1)
    III = SimpsonCompost(f, (a+b)/2, b, 1)
    Eab = abs(I - (II + III))
    if(Eab < eps*(b-a)):
        return I
    else:
        ABS.append((a+b)/2)
        return SimpsonAdaptativa2(f, a, (a+b)/2,eps, ABS) + SimpsonAdaptativa2(f, (a+b)/2, b, eps, ABS)

ABS1 = [a,b]
ABS2 = [a,b]
tol1 = 0.5e-3
tol2 = 0.5e-6
AP1 = SimpsonAdaptativa2(f, a, b, tol1, ABS1)
AP2 = SimpsonAdaptativa2(f, a, b, tol2, ABS2)
print("Usant un mètode recursiu i una tolerància de 1e-3 obtenim ", AP1)
print("Usant un mètode recursiu i una tolerància de 1e-6 obtenim ", AP2)
print("El primer resultat té un error absolut de: " , abs(AP1 - scipy.integrate.quad(f,a,b)[0]))
print("El segon resultat té un error absolut de: " , abs(AP2 - scipy.integrate.quad(f,a,b)[0]))
print(len(ABS2))
print(len(ABS1))

plt.title("Representació subintervals de la quadratura adaptativa amb la tol1")
plt.plot(np.unique(ABS1), f(np.unique(ABS1)), 'lightsteelblue')
plt.plot(np.unique(ABS1), f(np.unique(ABS1)), 'b.')
plt.legend()
plt.show()

plt.title("Representació subintervals de la quadratura adaptativa amb la tol2")
plt.plot(np.unique(ABS2), f(np.unique(ABS2)), 'springgreen')
plt.plot(np.unique(ABS2), f(np.unique(ABS2)), 'g.')
plt.show()
