function [ Smax, Smin ] = extremes_values( spot, rateclp, rateusd, vol, tenor,z )
Smax=spot*exp(((rateclp-rateusd)-((vol^2)/2))*tenor+vol*sqrt(tenor)*z);
Smin=spot*exp(((rateclp-rateusd)-((vol^2)/2))*tenor-vol*sqrt(tenor)*z);
end

