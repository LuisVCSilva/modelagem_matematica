from sympy import symbols,sympify,diff,integrate,Matrix,cos,sin	
from numpy import linspace
from matplotlib.pyplot import plot,show,xlabel,ylabel,legend

x,t,E,I,k,u_1,u_2,area_1,area_2,rho = \
symbols('x t E I k u_1 u_2 area_1 area_2 rho')

# polinomios de hermite de grau 3 e 4 - alternativamente chebysev
T_3 = sympify("4*x**3-3*x")
T_4 = sympify("8*x**4-8*x**2+1")  


#parametros de entrada
entrada_sistema = lambda t: cos(50*t)
comprimento_barra = 20

_k      = 1e+6
_rho    = 3.0
_area_1 = 0.9
_area_2 = 0.1


#matrizes de rigidez, inercia e massa
wxx_1 = Matrix([\
             [diff(diff(diff(T_3)))*T_3, diff(diff(diff(T_3)))*T_4],\
             [diff(diff(diff(T_4)))*T_3, diff(diff(diff(T_4)))*T_4]\
             ])

wxx_2 = Matrix([\
       [(diff(diff(T_3)))*diff(T_3), (diff(diff(T_3)))*diff(T_4)],\
       [(diff(diff(T_4)))*diff(T_3), (diff(diff(T_4)))*diff(T_4)]\
       ])

wxx_3 = Matrix([\
       [integrate(diff(diff(T_3))*diff(diff(T_3)),x),\
       integrate(diff(diff(T_4))*diff(diff(T_3)),x)],\
       [integrate(diff(diff(T_3))*diff(diff(T_4),x)),\
       integrate(diff(diff(T_4))*diff(diff(T_4)),x)]\
       ])


P1 = wxx_1-wxx_3.subs(x,comprimento_barra)-wxx_1-wxx_3.subs(x,0)

P2 = (-1/k)*wxx_2.subs(x,comprimento_barra)-wxx_2.subs(x,0)

P3 = rho*integrate(\
                   Matrix(\
                   [[T_3*T_3, T_3*T_4],\
                   [T_4*T_3, T_4*T_4]]\
                   ) *\
                   (area_1-area_2*x),(x,0,comprimento_barra))

P4 = -1*(entrada_sistema(t) * Matrix(\
                            [[T_3*T_3,T_3*T_4],\
                             [T_4*T_3,T_4*T_4]])\
        ).subs(x,comprimento_barra)-\
         (entrada_sistema(t)*Matrix([[T_3*T_3,T_3*T_4],\
                            [T_4*T_3,T_4*T_4]])\
         ).subs(x,0)

F = integrate(Matrix(\
              [[x*cos(t)*T_3],[x*cos(t)*T_4]]),\
              (x,0,1))

R1 = P4+P1

R2 = (P2+P3).subs((k, rho, area_1, area_2),(_k,_rho,_area_1,_area_2))
R2 = R2.subs(k, _k)
R2 = R2.subs(rho, _rho)
R2 = R2.subs(area_1, _area_1)
R2 = R2.subs(area_2, _area_2)

# Analise da dinamica de R1*u +R2*u''= F
barra_momento = R1*R2.inv()
barra_forca = R2.inv()*F


a = 0
b = 10
N = 10000
h = (b-a)/N

t = linspace(a,b,N,endpoint=True)

u1_v1 = [0.0] # deslocamento vertical
u1_v2 = [0.0] # velocidade vertical
u2_v1 = [0.0] # deslocamento horizontal
u2_v2 = [0.0] # velocidade horizontal

for i in range(len(t)-1):
    u1_v1.append(u1_v1[i]+h*u1_v2[i])
    u1_v2.append(eval(str(u1_v2[i]+h*((barra_momento[0,0].subs(t,t[i])*u1_v1[i]+\
              barra_momento[0,1].subs(t,t[i])*u2_v1[i])+\
              (F[0].subs(t,t[i]))).subs(t,t[i])\
                      ).replace("t",str(t[i])))\
             )
    u2_v1.append(u2_v1[i]+h*u2_v2[i])
    u2_v2.append(eval(str(u2_v2[i]+h*((barra_momento[1,0].subs(t,t[i])*u1_v1[i]+\
              barra_momento[1,1].subs(t,t[i])*u2_v1[i])+\
              (F[1].subs(t,t[i]))).subs(t,t[i])\
                      ).replace("t",str(t[i])))\
             )

plot(t,u1_v1)
plot(t,u2_v1)
plot(t,u1_v2)
plot(t,u2_v2)
legend(u1_v1,"$u_1$")
legend(u2_v1,"$u_2$")

legend(u1_v2,"${u'}_{1}$")
legend(u2_v2,"${u'}_{2}$")

xlabel("t")
ylabel("f(t)")
show()
