%% Electrolyte conductivity function: kappa(c_e) [1/Ohms*m]
%   Originally published by Scott Moura under GPL-3.0 license in December, 2016.
%   Link to original repository: https://github.com/scott-moura/SPMeT
%   Reorganized and extended by Nikolaos Wassiliadis, Andreas Bach, Matthias Wanzel in April, 2020.
%   Parameter taken from DUALFOIL LiPF6 in EC:DMC, Capiaglia et al. 1999

function [kappa] = electrolyteCond(c_e,T,p)

kappa = 0.0911+1.9101*c_e/1e3 - 1.052*(c_e/1e3).^2 + 0.1554*(c_e/1e3).^3;
kappa = kappa * exp(p.E.kappa_e/p.R*(1/p.T_ref_electrolyte - 1/T));
    
end