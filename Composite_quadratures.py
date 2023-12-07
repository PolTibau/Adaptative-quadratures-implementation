import numpy as np
from quadratures_simples import trobaPesos
from quadratures_simples import trobaPunts
from quadratures_Gauss import quadraturaGauss
 

def calcula(W,X,f):
    integral = 0
    for i in range(len(W)):
        for j in range(W[0].size):
            integral += W[i][j]*f(X[i][j])
    return integral


def trapeziCompost(f,a,b,m):
     h = (b-a)/m
     #Tindrem un vector de vectors de pesos per cada seccio W
     W = np.zeros((m,2))
     X = np.zeros((m,2))
     for i in range(m):
         X[i] = trobaPunts(1, a + i*h, a + (i+1)*h)
         W[i] = trobaPesos(1, a + i*h , a + (i+1)*h, X[i])
     return calcula(W,X,f)

def SimpsonCompost(f,a,b,m):
     h = (b-a)/m
     #Tindrem un vector de vectors de pesos per cada seccio W
     W = np.zeros((m,3))
     X = np.zeros((m,3))
     for i in range(m):
         X[i] = trobaPunts(2, a + i*h, a + (i+1)*h)
         W[i] = trobaPesos(2, a + i*h , a + (i+1)*h, X[i])
     return calcula(W,X,f)

def GaussCompost(f,a,b,m,nIP):
    h = (b-a)/m
    W = np.zeros(nIP,m)
    for i in range(m):
        W[i] = quadraturaGauss(nIP)
        X[i] = trobaPunts(1, a + i*h, a + (i+1)*h)

