function [value_ddmc_omega,accuracy_omega] = DDMC_omega( spot, rateclp, rateusd, spotbarra, strike, tenor, N, dt, eta )
rng default %semilla
z=randn(N,tenor);%Se genera un numero random de tamano N, tenor
value_ddmc_omega=zeros(1,21);
accuracy_omega=zeros(1,21);
k=0;
for omega=0:0.05:1
    k=k+1;
    S(1)=spot; %El primer valor que tomara S sera 'spot'
        for j=1:N %Genera de 1 a N simulaciones
            for i=1:tenor
                S(i+1)=S(i)*exp(((rateclp-rateusd-((eta*(omega+(1-omega).*(spotbarra./S(i))))^2)/2)).*dt+(eta*(omega+(1-omega).*(spotbarra./S(i)))).*sqrt(dt).*z(j,i));
            end
            VP(j)=exp(-rateclp.*tenor)*max(S(end)-strike,0); %Payoffs descontados a valor presente
        end
a0=N;       %Calculo de acumuladores
a1=sum(VP);
a2=sum(VP.^2);
value_ddmc_omega(1,k)=(a1/a0);
accuracy_omega(1,k)=(1/sqrt(a0))*sqrt((a2/a0)-(value_ddmc_omega(1,k)).^2);
end
end
