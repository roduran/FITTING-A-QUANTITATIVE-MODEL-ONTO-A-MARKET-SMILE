function [ A_DD ] = MatrizA_DD( spot, rateclp, rateusd, tenor, N, vol, z, eta, omega, spotbarra)
A_DD=zeros(N,N);
[ Smax, Smin ] = extremes_values( spot, rateclp, rateusd, vol, tenor,z );
S=linspace(Smin,Smax);
dS=S(2)-S(1);
d(1)=((rateclp-rateusd)*S(1))/dS;
d(N)=((rateclp-rateusd)*S(N))/dS;

A_DD(1,1)=-d(1)-rateclp;
A_DD(1,2)=d(1);
A_DD(N,N-1)=-d(N);
A_DD(N,N)=d(N)-rateclp;

for j=2:N-1
    d(j)=((rateclp-rateusd)*S(j))/dS;
    g(j)=((eta*(omega+(1-omega).*(spotbarra./S(j))))^2)*(S(j)^2)/dS^2;
    
    A_DD(j,j-1)=-d(j)/2+g(j)/2;
    A_DD(j,j)=-g(j)-rateclp;
    A_DD(j,j+1)=d(j)/2+g(j)/2;
end
end


