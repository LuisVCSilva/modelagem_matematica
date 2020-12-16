from matplotlib.pyplot import *
from numpy import *

def FiltroGaussiano (t, freq_central=5, larg_banda=0.5, \
                     ref_larg_banda=-6,limiar_truncamento=-60):
   #frequencia central deve ser maior ou igual a zero
   #largura de banda fracionaria deve ser maior que zero
   #limiar de referencia para largura de banda deve ser menor que 0
   #ja que a funcao gaussiana eh para x \in (-\infty,\infty), o limiar de 
   #truncamento trunca os valores do filtro apos um ponto

   # exp(-a t^2) <->  sqrt(pi/a) exp(-pi^2/a * f^2)  = g(f)
   ref = pow(10.0, ref_larg_banda / 20.0)
   # fdel = freq_central*larg_banda/2:  g(fdel) = ref
   # colocando em funcao de a temos que...
   # pi^2/a * freq_central^2 * larg_banda^2 /4=-log(ref)
   a = -(pi * freq_central * larg_banda) ** 2 / (4.0 * log(ref))

   sinal_envoltorio = exp(-a * t * t)
   sinal_real = sinal_envoltorio * cos(2 * pi * freq_central * t)
   sinal_imaginario = sinal_envoltorio * sin(2 * pi * freq_central * t)
   return sinal_real, sinal_imaginario, sinal_envoltorio


intervalo = (-1.0,1.0)
h = 0.001

t = linspace(intervalo[0], intervalo[1], 1/h)
freq_central       =  4.0
larg_banda         =  0.5
ref_larg_banda     = -6.0
limiar_truncamento = -60.0

sinal_real, sinal_imaginario, envoltorio = FiltroGaussiano(
                                             t, \
                                             freq_central, \
                                             larg_banda, \
                                             ref_larg_banda, \
                                             limiar_truncamento \
                                             )
plot(t, sinal_real, t, sinal_imaginario, t, envoltorio, '--')

legend(["Parte Real","Parte Imaginaria","Envelope"])
title("Constante de Tempo RC de Filtro Gaussiano\n" + \
      "$e^{-at^2}e^{1j 2\pi f_c t }$")
show()
