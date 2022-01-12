%% Cathode reaction rate function: k_p(SOC_p) [(A/m^2)*(mol^3/mol)^(1+alpha)]
%   Called from ode_spmet.m  
%   Originally published by Scott Moura under GPL-3.0 license in December, 2016.
%   Link to original repository: https://github.com/scott-moura/SPMeT
%   Extended with a SOC- and T-dependent diffusion and reaction rate for both electrodes by Nikolaos
%   Wassiliadis, Andreas Bach, Matthias Wanzel in April, 2020.

function [k_p] = cathodeK(SOC_p, p)
        k_p = p.kin_k_p(1-SOC_p) * p.K_multi_p;
end
