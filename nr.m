function [ sigmanew ] = nr( vol, ValorObjetivo, ValorBS, VegaBS )

sigmanew=vol + (ValorObjetivo - ValorBS)/VegaBS;

end

