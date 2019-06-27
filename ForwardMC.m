function [value_f,accuracy ] = ForwardMC( spot, rateclp, rateusd, vol, strike, tenor, N,dt )
rng default %semilla
S(1)=spot; %El primer valor que tomara S sera 'spot'
z=randn(N,1/dt);%Se genera un numero random de tamano N, tenor
for j=1:N %Genera de 1 a N simulaciones
    for i=1:1/dt
        S(i+1)=S(i)*exp(((rateclp-rateusd-(vol^2)/2)).*dt+vol.*sqrt(dt).*z(j,i));
    end
M=max(S(end)-strike,0)-max(strike-S(end),0); %Calculo valor Straddle
VP(j)=exp(-rateclp.*tenor).*M; %Payoffs descontados a valor presente
end

a0=N;       %Calculo de acumuladores
a1=sum(VP);
a2=sum(VP.^2);

value_f=a1/a0;
accuracy=(1/sqrt(a0))*sqrt((a2/a0)-(value_f).^2);

end


