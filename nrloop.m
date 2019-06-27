function [ y ] = nrloop( vol, strike, tenor, spot, rateclp, rateusd, e, N, ValorObjetivo )

y=zeros(1,5);

for j=1:5
    sigma=vol;
    for i= 1:N
        [ ValorBS, VegaBS ] = Pr1(strike(1,j), tenor, spot, rateclp, rateusd, sigma, e);
        [ sigma ] = nr( sigma, ValorObjetivo(1,j), ValorBS, VegaBS );
    end

    y(1,j)= sigma;
end
end

