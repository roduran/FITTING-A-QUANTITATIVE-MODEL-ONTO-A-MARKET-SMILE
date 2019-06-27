function [ ValorBS, VegaBS ] = Pr1(strike, tenor, spot, rateclp, rateusd, sigma, e)

d1=((log(spot/strike)+(rateclp-rateusd)*tenor)/(sigma*sqrt(tenor))) + (sigma*sqrt(tenor)/2); 
d2=d1-sigma*sqrt(tenor);

ValorBS= e*spot*exp(-rateusd*tenor)*normcdf(e*d1) - e*strike*exp(-rateclp*tenor)*normcdf(e*d2); 

VegaBS= spot * normpdf(d1) * sqrt(tenor)*exp(-rateusd*tenor);
end
