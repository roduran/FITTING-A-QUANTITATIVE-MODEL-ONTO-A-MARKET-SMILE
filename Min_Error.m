function [Error] = Min_Error(x,i)

% Data
load('K.mat');
load('M1.mat');
load('M2.mat');
load('M3.mat');
load('M6.mat');
load('Y1.mat');
load('S.mat');
load('domestic.mat');
load('foreign.mat');
load('T.mat');
% Actualizacion de rates
for b=1:5
    for p = 1:2663
        rateclp(p,b) = -log(domestic(p,b))/(T(p,b)/365);
        rateusd(p,b) = -log(foreign(p,b))/(T(p,b)/365);
    end
end
% Parametros
VolatilidadMercado = [M1 M2 M3 M6 Y1]/100;
z =3;
N=100;
dt =1/365;
% Entrega posicion parametros iniciales
omega0 = x(1);
eta0 = x(2);
% Fraccionamiento de la Data
for j=1:5
     if j==1
        h=1;
        H=5;
    end
    if j==2
        h=6;
        H=10;
    end
    if j==3
        h=11;
        H=15;
    end
    if j==4
        h=16;
        H=20;
    end
    if j==5
        h=21;
        H=25;
    end

for u=h:H
    
    [vol0] = DisplacedDiffusion(eta0,omega0,S(i));
    [ value_cn ] = DF_CN_DD(S(i),K(i,u),rateclp(i,j),rateusd(i,j),vol0,T(i,j),z,N,dt,eta0,omega0);
    Volatilidad(1,u) = NewtonRaphson(K(i,u),S(i),rateclp(i,j),rateusd(i,j),T(i,j),value_cn, VolatilidadMercado(i,u));
    Vol_Implicita(1,u)=Volatilidad(u);
    
end

end
error = abs(VolatilidadMercado(i,:)-Vol_Implicita);
Error = mean(error);
end