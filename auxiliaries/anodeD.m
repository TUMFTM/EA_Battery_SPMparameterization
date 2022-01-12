%% Anode diffusion coefficient function: D_S_n(SOC_n) [m^2/s]
%   Called by ode_spmet.m
%   Originally published by Scott Moura under GPL-3.0 license in December, 2016.
%   Link to original repository: https://github.com/scott-moura/SPMeT
%   Extended with a SOC- and T-dependent diffusion and reaction rate for both electrodes by Nikolaos
%   Wassiliadis, Andreas Bach, Matthias Wanzel in April, 2020.

function [D_s_n] = anodeD(SOC_n,p)
       D_s_n = p.kin_D_n(SOC_n) * p.D_multi_n;
end
