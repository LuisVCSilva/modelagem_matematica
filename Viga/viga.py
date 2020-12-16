from sympy import symbols,sympify,diff,integrate,Matrix,cos,sin
from numpy import linspace
from matplotlib.pyplot import plot,show,xlabel,ylabel

x,t,E,I,k,u_1,u_2,delta_1,delta_2,rho = symbols('x t E I k u_1 u_2 delta_1 delta_2 rho')

t3 = sympify("4*x**3-3*x")       # polinomios de hermite de grau 3 e 4 - alternativamente chebysev
t4 = sympify("8*x**4-8*x**2+1")  

wxx_1 = Matrix([\
             [diff(diff(diff(t3)))*t3,\
             diff(diff(diff(t3)))*t4],\
             [diff(diff(diff(t4)))*t3,\
             diff(diff(diff(t4)))*t4]\
      ])

wxx_2 = Matrix([\
       [(diff(diff(t3)))*diff(t3),\
       (diff(diff(t3)))*diff(t4)],\
       [(diff(diff(t4)))*diff(t3),\
       (diff(diff(t4)))*diff(t4)]\
      ])

wxx_3 = Matrix([\
       [integrate(diff(diff(t3))*diff(diff(t3)),x),\
       integrate(diff(diff(t4))*diff(diff(t3)),x)],\
       [integrate(diff(diff(t3))*diff(diff(t4),x)),\
       integrate(diff(diff(t4))*diff(diff(t4)),x)]\
      ])

intervalo_integracao = 15
P1 = wxx_1-wxx_3.subs(x,intervalo_integracao)-wxx_1-wxx_3.subs(x,0)

P2 = (-1/k)*wxx_2.subs(x,intervalo_integracao)-wxx_2.subs(x,0)

P3 = rho*integrate(Matrix([[t3*t3, t3*t4], [t4*t3, t4*t4]]) * (delta_1-delta_2*x),(x,0,intervalo_integracao))

entrada_sistema = lambda t: cos(50*t)

P4 = -1*(entrada_sistema(t) * Matrix([[t3*t3,t3*t4],[t4*t3,t4*t4]])).subs(x,intervalo_integracao)\
        -(entrada_sistema(t)*Matrix([[t3*t3,t3*t4],[t4*t3,t4*t4]])).subs(x,0)

F = integrate(Matrix([[x*cos(t)*t3],[x*cos(t)*t4]]),(x,0,1))

R1 = P4+P1

R2 = (P2+P3).subs((k, rho, delta_1, delta_2),(10**6, 3, 0.5, 0.2))
R2 = R2.subs(k, 10**6)
R2 = R2.subs(rho, 3)
R2 = R2.subs(delta_1, 0.5)
R2 = R2.subs(delta_2, 0.2)

# Analise da dinamica de R1*u +R2*u'' = F
barra_momento = R1*R2.inv()
barra_forca = R2.inv()*F


a = 0
b = 10
N = 1000
h = (b-a)/N

t = linspace(a,b,N,endpoint=True)

v1 = [0.0]
v2 = [0.0]
v3 = [0.0]
v4 = [0.0]

for i in range(len(t)-1):
    v1.append(v1[i]+h*v2[i])
    v2.append(eval(str(v2[i]+h*((barra_momento[0,0].subs(t,t[i])*v1[i]+barra_momento[0,1].subs(t,t[i])*v3[i])+(F[0].subs(t,t[i]))).subs(t,t[i])).replace("t",str(t[i]))))
    v3.append(v3[i]+h*v4[i])
    v4.append(eval(str(v4[i]+h*((barra_momento[1,0].subs(t,t[i])*v1[i]+barra_momento[1,1].subs(t,t[i])*v3[i])+(F[1].subs(t,t[i]))).subs(t,t[i])).replace("t",str(t[i]))))

plot(t,v1)
plot(t,v3)
plot(t,v2)
plot(t,v4)
xlabel("t")
ylabel("f(t)")
show()
