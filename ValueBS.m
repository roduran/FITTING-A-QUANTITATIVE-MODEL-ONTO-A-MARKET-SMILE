function[valueBS,vega] = ValueBS(K,S,r,q,T,vol0)
    e = 1;

    t=T/365;
    d1 = (log(S/K)+(r-q)*t)/(vol0*sqrt(t)) + vol0*sqrt(t)/2;
    d2 = d1 - vol0*sqrt(t);

    valueBS = e*S*exp(-q*t)*normcdf(e*d1)-e*K*exp(-r*t)*normcdf(e*d2);
    vega = S*normpdf(d1)*sqrt(t);

end