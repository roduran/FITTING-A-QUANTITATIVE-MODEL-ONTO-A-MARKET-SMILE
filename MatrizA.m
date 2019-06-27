function [ A ] = MatrizA( spot, rateclp, rateusd, tenor, N, vol, z)
A=zeros(N,N);
[ Smax, Smin ] = extremes_values( spot, rateclp, rateusd, vol, tenor,z );
S=linspace(Smin,Smax);
dS=S(2)-S(1);
d(1)=((rateclp-rateusd)*S(1))/dS;
d(N)=((rateclp-rateusd)*S(N))/dS;

A(1,1)=-d(1)-rateclp;
A(1,2)=d(1);
A(N,N-1)=-d(N);
A(N,N)=d(N)-rateclp;

for j=2:N-1
    d(j)=((rateclp-rateusd)*S(j))/dS;
    g(j)=(vol^2)*(S(j)^2)/dS^2;
    
    A(j,j-1)=-d(j)/2+g(j)/2;
    A(j,j)=-g(j)-rateclp;
    A(j,j+1)=d(j)/2+g(j)/2;
end
end

