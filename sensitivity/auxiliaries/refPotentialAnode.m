%% Reference potential for anode: Unref(theta_n)
%   Called by ode_spmet.m
%   Originally published by Scott Moura under GPL-3.0 license in December, 2016.
%   Link to original repository: https://github.com/scott-moura/SPMeT
%   Reorganized and extended by Nikolaos Wassiliadis, Andreas Bach, Matthias Wanzel in April, 2020.

function [Uref,varargout] = refPotentialAnode(p,theta,pscale)

theta = real(theta);
if theta > 1
    theta = 1;
elseif theta < 0
    theta = 0 ;
end

if(~isreal(theta))
    beep;
    error('dfn:err','Complex theta_n');
end

Uref = (p.p_n(1)+p.p_n(2)*exp(p.p_n(3)*theta)+p.p_n(4)*tanh((theta+p.p_n(5))/p.p_n(6)) ... 
       +p.p_n(7)*tanh((theta+p.p_n(8))/p.p_n(9))+p.p_n(10)*tanh((theta+p.p_n(11))/p.p_n(12)) ... 
       +p.p_n(13)*tanh((theta+p.p_n(14))/p.p_n(15))); 
 
% Gradient of OCP wrt theta
if(nargout >= 2)
    varargout{1} = 0;
end

% Gradient of OCP wrt temperature
if(nargout >= 3)   
    varargout{2} = 0;   
end

