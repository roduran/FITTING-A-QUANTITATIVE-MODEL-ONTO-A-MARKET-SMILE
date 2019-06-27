function [value_cn] = DF_CN_DD(S,K,rateclp,rateusd,vol0,T,z,N,dt,eta0,omega0)
%SpaceNodes
t=T/365;
Smin = S*exp((rateclp-rateusd-(vol0^2)/2).*(t)-vol0*sqrt(t).*z);
Smax = S*exp((rateclp-rateusd-(vol0^2)/2).*(t)-vol0*sqrt(t).*(-z));
Spot = linspace(Smin,Smax,N);
V = max(Spot-K,0);
Vt = V';
%Matriz A
A = zeros(N,N);
dS = Spot(2)-Spot(1);
d(1) = ((rateclp-rateusd)*Spot(1))/dS;
d(100) = ((rateclp-rateusd)*Spot(100))/dS;
A(1,1) = -d(1)-rateclp;
A(1,2) = d(1);
A(100,99) = -d(100);
A(100,100) = d(100)-rateclp;

for i = 2:N-1 
    d(i) = ((rateclp-rateusd)*Spot(i))/dS;
    %Modelo VL
    g(i) = ((eta0*(omega0+(1-omega0).*(600/Spot(i))))^2 * Spot(i)^2)/dS^2;
    A(i,i) = -g(i)-rateclp;
    A(i,i-1) = -d(i)/2 + g(i)/2;
    A(i,i+1) = d(i)/2 + g(i)/2; 
end
%Valorizacion
I = eye(N);
R=inv(I-A*dt/2);
RR=I+A*dt/2;
    
for i = 1:t/dt 
    Vt = (R)*(RR)*Vt;    
end 

value_cn = interp1(Spot,Vt,S);
end