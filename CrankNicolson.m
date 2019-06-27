function [ value_cn, vector_cn ] = CrankNicolson( spot, strike, rateclp, rateusd, vol, tenor, z, N, dtau )
I=eye(N);
vector_cn=zeros(N,1);
[ A ] = MatrizA( spot, rateclp, rateusd, tenor, N, vol, z);
[ Smax, Smin ] = extremes_values( spot, rateclp, rateusd, vol, tenor,z );
[ S ] = SpaceNodes( Smax, Smin, N );
value_call= (max(S-strike,0))';

for i=dtau:dtau:tenor
    value_call=(inv(I-A*dtau/2))*(I+A*dtau/2)*value_call;
end

value_cn=interp1(S, value_call, spot);

end


