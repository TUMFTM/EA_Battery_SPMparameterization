%% Cathode diffusion coefficient function: D_s_p(SOC_p) [m^2/s]
%   Called by ode_spmet.m  
%   Originally published by Scott Moura under GPL-3.0 license in December, 2016.
%   Link to original repository: https://github.com/scott-moura/SPMeT
%   Extended with a SOC- and T-dependent diffusion and reaction rate for both electrodes by Nikolaos
%   Wassiliadis, Andreas Bach, Matthias Wanzel in April, 2020.

function [D_s_p] = cathodeD(SOC_p, p)
        D_s_p = p.kin_D_p(1-SOC_p) * p.D_multi_p;
end
