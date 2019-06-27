function [Volatilidad] = NewtonRaphson(K,S,rateclp,rateusd,T,value_cn, VolatilidadMercado)
    accuracy = 10.^(-2);
    acc = 1+accuracy;
    Volatilidad = 0.1;

    while abs(acc)>accuracy 
        vol0 = Volatilidad;
        [valueBS,vega] = ValueBS(K,S,rateclp,rateusd,T,vol0);
        Volatilidad = vol0 + (value_cn-valueBS)/vega;
        acc = abs(Volatilidad-vol0);
    end
     if isnan(Volatilidad)
      Volatilidad=VolatilidadMercado;  
    end
    end
    

