function [ value ] = EuropeanCallOption( Spot, r, q, vol, Strike, T )

d1=(log(Spot/Strike)+(r-q)*T)/(vol*sqrt(T))+(vol*sqrt(T))/2;
d2=d1-vol*sqrt(T);
N1=normcdf(d1);
N2=normcdf(d2);

value=(Spot*exp(-q*T)*N1)-(Strike*exp(-r*T)*N2);

end

