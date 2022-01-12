%% Reference potential for cathode: Upref(theta_p)
%   Called from ode_spmet.m
%   Originally published by Scott Moura under GPL-3.0 license in December, 2016.
%   Link to original repository: https://github.com/scott-moura/SPMeT
%   Reorganized and extended by Nikolaos Wassiliadis, Andreas Bach, Matthias Wanzel in April, 2020.

function [Uref,varargout] = refPotentialCathode(p,theta,pscale)

theta = real(theta);
if theta > 1
    theta = 1;
elseif theta < 0
    theta = 0 ;
end

Uref =  (p.p_p(1)*(1-theta).^14   + p.p_p(2)*(1-theta).^13 + p.p_p(3)*(1-theta).^12 ... 
        + p.p_p(4)*(1-theta).^11 + p.p_p(5)*(1-theta).^10 + p.p_p(6)*(1-theta).^9 ... 
        + p.p_p(7)*(1-theta).^8  + p.p_p(8)*(1-theta).^7  + p.p_p(9)*(1-theta).^6 ...
        + p.p_p(10)*(1-theta).^5 + p.p_p(11)*(1-theta).^4 + p.p_p(12)*(1-theta).^3 ... 
        + p.p_p(13)*(1-theta).^2 + p.p_p(14).*(1-theta)   + p.p_p(15));
    
% Gradient of OCP wrt theta
if(nargout == 2)
    varargout{1} = 0;   
end

% Gradient of OCP wrt temperature
if(nargout >= 3)  
    varargout{2} = 0;   
end