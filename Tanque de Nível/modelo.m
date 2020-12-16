clc, clear all, close all

A = 10;            % area do tanque em m^2
H(1) = 30;         % altura inicial da coluna de liquido no tanque
a = 0;             % ponto inicial usado no metodo de euler
b = 5;             % ponto final usado no metodo de euler
N = 100;           % numero de pontos usado no metodo de euler
passo = (b-a) / N; % passo usado no método de euler
t = passo.*(0:1:N) % tempo de 0 a 5 segundos
Qi = 20.*t;        % vazão de entrada m^3/s
Qo = 10.*t;        % vazão de saída  m^3/s
 
for i = 1:N
    H(i+1) = H(i) + passo.*(Qi(i) - Qo(i))./ A ; % altura em função dos parametros de entrada
end
 
plot(t,H);
xlabel('Tempo (s)');
ylabel('Nivel de fluido no tanque(m)');
grid();
sleep(3);