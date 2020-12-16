from matplotlib.pyplot import *
from numpy import *

A = 10.0             # area do tanque (m**2)
H = [100.0]           # altura inicial de fluido no tanque (m)
a = 0.0              # intervalos de integracao
b = 10.0              # 
N = 1000              # discretizacao do intervalo de integracao
passo = (b-a) / N    # passo
t = linspace(\
             a,\
             b,\
             N,\
             endpoint=False\
            ) # vetor de tempo

fQe = lambda t: 10*t # razao da vazao de entrada (m**3)
fQs = lambda t: 50*pi*cos(t) # razao de vazao de saida   (m**3)

Qe = fQe(t)         
Qs = fQs(t)

for i in range(N-1):
    H.append(H[i] + passo*(Qe[i] - Qs[i])/A)  # altura em função dos parametros de entrada

 
plot(t.tolist(),H)
xlabel('Tempo (s)')
ylabel('Altura de fluido no tanque (m)')
grid()
show()

