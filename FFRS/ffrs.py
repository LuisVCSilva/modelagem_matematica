from __future__ import print_function
from itertools import product as produto_cartesiano
from matplotlib.pyplot import *
from scipy import signal
from random import randint as rand
from numpy import *

def FFRS(r, s,q):#considera-se -1 sinal invalido
     return 0 if (r,s,q) in [(0,0,0),(0,1,0),(0,1,1)] else (1 if (r,s,q) in [(0,0,1),(1,0,0),(1,0,1)] else -1)



def tabela_verdade():
 print('r\ts\tQ(k)\tQ(k+1)')
 for r, s,q in produto_cartesiano((0, 1), repeat=3):
    print('{!s:5}\t{!s:5}\t{!s:5}\t{!s:5}'.format(r, s,q, FFRS(r,s,q)))

print("Tabela Verdade do Flip Flop RS\n"+\
      "Considere \"-1\" como estado proibido\n"\
     )
tabela_verdade()

tamanho_palavra = 1024

t = linspace(0, 1, tamanho_palavra, endpoint=True)

r = [0 if x<0 else 1 for x in signal.square(5 * pi * t)]
subplot(4, 1, 1)
plot(t, r)

s = [0 if x<0 else 1 for x in signal.square(10 * pi * t)]
subplot(4, 1, 2)
plot(t, s)

q = [0 if x<0 else 1 for x in signal.square(15 * pi * t)]
subplot(4, 1, 3)
plot(t, q)

fx = [FFRS(_r,_s,_q) for (_r,_s,_q) in zip(r,s,q)]
subplot(4, 1, 4)
plot(t, fx)

xlabel("t")
ylabel("f(t)")

ylim(-1.5, 1.5)
plot()
show()
