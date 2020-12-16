clc;clear all;close all;

syms x t E I k u1 u2 delta1 delta2 rho

t3 = 4*x^3-3*x;
t4 = 8*x^4-8*x^2+1;

% ppp1=[diff(diff(diff(t3)))*t3 diff(diff(diff(t3)))*t4;diff(diff(diff(t4)))*t3 diff(diff(diff(t4)))*t4];
% ppp2=[(diff(diff(t3)))*diff(t3) (diff(diff(t3)))*diff(t4);(diff(diff(t4)))*diff(t3) (diff(diff(t4)))*diff(t4)];
% ppp3=[(diff(diff(t3)))*diff(t3) (diff(diff(t4)))*diff(t3);(diff(diff(t3)))*diff(t4) (diff(diff(t4)))*diff(t4)];
% ppp4=[diff(diff(diff(t3)))*t3 diff(diff(diff(t4)))*t3;diff(diff(diff(t3)))*t4 diff(diff(diff(t4)))*t4];
% 
% P1=subs(vpa(E*I*(ppp1-ppp2-ppp3+ppp4)),x,1)-subs(vpa(E*I*(ppp1-ppp2-ppp3+ppp4)),x,0);

pp1 = [diff(diff(diff(t3)))*t3 diff(diff(diff(t3)))*t4;diff(diff(diff(t4)))*t3 diff(diff(diff(t4)))*t4];
pp2 = [(diff(diff(t3)))*diff(t3) (diff(diff(t3)))*diff(t4);(diff(diff(t4)))*diff(t3) (diff(diff(t4)))*diff(t4)];
pp3 = [int(diff(diff(t3))*diff(diff(t3))) int(diff(diff(t4))*diff(diff(t3)));int(diff(diff(t3))*diff(diff(t4))) int(diff(diff(t4))*diff(diff(t4)))];

l = 15;

P1 = subs(pp1-pp3,x,l)-subs(pp1-pp3,x,0);

P2 = (-1/k)*subs(pp2,x,l)-subs(pp2,x,0);

P3 = rho*int([t3*t3 t3*t4;t4*t3 t4*t4]*(delta1 -delta2*x),0,l);

P4 = -1*(subs(cos(20*t)*[t3*t3 t3*t4;t4*t3 t4*t4],x,l)-subs(cos(20*t)*[t3*t3 t3*t4;t4*t3 t4*t4],x,0));

F = int([x*cos(t)*t3;x*cos(t)*t4],0,1);

R1 = P4+double(P1);

R2 = subs(P2+P3,{k rho delta1 delta2},[10^6 3 0.5 0.2]);

%R1*u+R2*ü=F

MBARRA = R1*double(inv(double(R2)));
FBARRA = inv(double(R2))*F;


a=0;
b=10;
n=100;
h=(b-a)/n;

t=[a:h:b];

v1(1)=0;
v2(1)=0;
v3(1)=0;
v4(1)=0;

for i=1:length(t)-1
    v1(i+1)=v1(i)+h*v2(i);
    v2(i+1)=v2(i)+h*(subs(MBARRA(1,1),t(i))*v1(i)+subs(MBARRA(1,2),t(i))*v3(i))+subs(F(1),t(i));
    v3(i+1)=v3(i)+h*v4(i);
    v4(i+1)=v4(i)+h*(subs(MBARRA(2,1),t(i))*v1(i)+subs(MBARRA(2,2),t(i))*v3(i))+subs(F(2),t(i));
end

plot(t,v1);
hold on
plot(t,v3);
figure,
plot(t,v2);
hold on
plot(t,v4);