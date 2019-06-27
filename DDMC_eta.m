function [value_ddmc_eta,accuracy ] = DDMC_eta( spot, rateclp, rateusd, spotbarra, strike, tenor, N, dt, omega )
rng default %semilla
S(1)=spot; %El primer valor que tomara S sera 'spot'
z=randn(N,tenor);%Se genera un numero random de tamano N, tenor
for eta=1:11
for j=1:N %Genera de 1 a N simulaciones
    for i=1:tenor
        S(i+1)=S(i)*exp(((rateclp-rateusd-((eta*(omega+(1-omega).*(spotbarra./spot)))^2)/2)).*dt+(eta*(omega+(1-omega).*(spotbarra./spot))).*sqrt(dt).*z(j,i));
    end
M=max(S(end)-strike,0); %Calculo valor Call
VP(j)=exp(-rateclp.*tenor).*M; %Payoffs descontados a valor presente
end

a0=N;       %Calculo de acumuladores
a1=sum(VP);
a2=sum(VP.^2);

value_ddmc_eta=a1/a0;
accuracy=(1/sqrt(a0))*sqrt((a2/a0)-(value_ddmc_eta).^2);
end
end

