%% Anode reaction rate function: k_n(SOC_n) [(A/m^2)*(mol^3/mol)^(1+alpha)] 
%   Called by ode_spmet.m 
%   Originally published by Scott Moura under GPL-3.0 license in December, 2016.
%   Link to original repository: https://github.com/scott-moura/SPMeT
%   Extended with a SOC- and T-dependent diffusion and reaction rate for both electrodes by Nikolaos
%   Wassiliadis, Andreas Bach, Matthias Wanzel in April, 2020.

function [k_n] = anodeK(SOC_n, p)
        k_n = p.kin_k_n(SOC_n) * p.K_multi_n;
end
