function [ value_f_cn ] = Forward_CN( spot, strike, rateclp, rateusd, vol, tenor, z, N, dtau )
I=eye(N);
[ A ] = MatrizA( spot, rateclp, rateusd, tenor, N, vol, z);
[ Smax, Smin ] = extremes_values( spot, rateclp, rateusd, vol, tenor,z );
[ S ] = SpaceNodes( Smax, Smin, N );
value_call= (max(S-strike,0)-max(strike-S,0))';

for i=dtau:dtau:tenor
    value_call=(inv(I-A*dtau/2))*(I+A*dtau/2)*value_call;
end

value_f_cn=interp1(S, value_call, spot);

end

