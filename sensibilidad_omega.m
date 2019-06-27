function [ value_DFDD_omega ] = sensibilidad_omega( spot, strike, rateclp, rateusd, vol, tenor, z, N, dtau, eta, spotbarra )
k=0;
value_DFDD_omega=zeros(1,21);
for omega=0:0.05:1
    k=k+1;
    [ A_DD ] = MatrizA_DD( spot, rateclp, rateusd, tenor, N, vol, z, eta, omega, spotbarra);
    [ value_cn_DD ] = CrankNicolson_DD( spot, strike, rateclp, rateusd, vol, tenor, z, N, dtau, eta, omega, spotbarra );
    value_DFDD_omega(1,k)=value_cn_DD;
end
end

