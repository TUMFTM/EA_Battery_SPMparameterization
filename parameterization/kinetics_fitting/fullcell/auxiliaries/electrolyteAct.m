%% Electrolyte activity coefficient function: dlnf/dln(c_e) 
%   Called by ode_spmet.m
%   Originally published by Scott Moura under GPL-3.0 license in December, 2016.
%   Link to original repository: https://github.com/scott-moura/SPMeT
%   Reorganized and extended by Nikolaos Wassiliadis, Andreas Bach, Matthias Wanzel in April, 2020.
%   Regression taken from Lundgren et al. (2014): DOI:10.1149/2.0641503jes.

function [dActivity] = electrolyteAct(c_e)

c_e = c_e./1000; %UnitConversion: 1 mol/L -> 1000 mol/m^3

dActivity =  (0.2731.*c_e.^2+0.6352.*c_e+0.4577)./(0.1291.*c_e.^3-0.3517.*c_e.^2+0.4893.*c_e+0.5713)-1;

end