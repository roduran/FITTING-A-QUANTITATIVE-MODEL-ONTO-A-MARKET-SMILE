function [vol0] = DisplacedDiffusion(eta0,omega0,S)
vol0 = eta0*(omega0+(1-omega0)*(600/S));
end

