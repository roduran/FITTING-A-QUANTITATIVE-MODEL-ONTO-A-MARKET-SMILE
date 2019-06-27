function [ value_cn_DD ] = CrankNicolson_DD( spot, strike, rateclp, rateusd, vol, tenor, z, N, dtau, eta, omega, spotbarra )
I=eye(N);
[ A_DD ] = MatrizA_DD( spot, rateclp, rateusd, tenor, N, vol, z, eta, omega, spotbarra);
[ Smax, Smin ] = extremes_values( spot, rateclp, rateusd, vol, tenor,z );
[ S ] = SpaceNodes( Smax, Smin, N );
value_call= (max(S-strike,0))';

for i=dtau:dtau:tenor
    value_call=(inv(I-A_DD*dtau/2))*(I+A_DD*dtau/2)*value_call;
end

value_cn_DD=interp1(S, value_call, spot);

end


