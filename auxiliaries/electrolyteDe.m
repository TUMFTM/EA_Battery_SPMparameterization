%% Electrolyte diffusion coefficient function: D_e(c_e,T) [m^2/s]
%   Called from ode_spmet.m
%   Originally published by Scott Moura under GPL-3.0 license in December, 2016.
%   Link to original repository: https://github.com/scott-moura/SPMeT
%   Reorganized and extended by Nikolaos Wassiliadis, Andreas Bach, Matthias Wanzel in April, 2020.
%   Parameter taken from Albertus et al. (2009) DOI:10.1149/1.3129656 and Ecker et al. (2015) DOI:
%   10.1149/2.0541509jes.


function [D_e,varargout] = electrolyteDe(c_e)

D_e = 6.5e-10*exp(-0.7*c_e/1e3);

if(nargout == 2)
    
    dD_e = -0.7*D_e/1e3;
    varargout{1} = dD_e;
end
end