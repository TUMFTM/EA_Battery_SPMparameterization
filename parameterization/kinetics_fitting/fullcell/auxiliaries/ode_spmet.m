%% ODEs
%   Called by spmet.m
%   Originally published by Scott Moura under GPL-3.0 license in December, 2016.
%   Link to original repository: https://github.com/scott-moura/SPMeT
%   Restructrued and extended by Nikolaos Wassiliadis, Andreas Bach, Matthias Wanzel in April, 2020.

function [x_dot,varargout] = ode_spmet(t,x,data,p)


%% Parse Input Data

% Parse and interpolate current
cur = griddedInterpolant(data.time,data.cur);
cur = cur(t);

% Parse states
c_s_n = real(x(1:(p.Nr-1)));
c_s_p = real(x(p.Nr : 2*(p.Nr-1)));
c_e = real(x(2*p.Nr-1 : 2*p.Nr-1+p.Nx-4));
delta_sei = 0;
T1 = p.T_amb; % Iso-thermal

% Calculate SOC
r_vec = (0:p.delta_r_n:1)';
c_n = [c_s_n(1); c_s_n; c_s_n(p.Nr-1)];
c_p = [c_s_p(1); c_s_p; c_s_p(p.Nr-1)];
SOC_n = 3/p.c_s_n_max * trapz(r_vec,r_vec.^2.*c_n);
SOC_p = 3/p.c_s_p_max * trapz(r_vec,r_vec.^2.*c_p);

%% Pre-calculations with current states

%%% MOLAR FLUXES
% Compute total molar flux
jn_tot = -cur/(p.Faraday*p.a_s_n*p.Area*p.L_n);
jp_tot = cur/(p.Faraday*p.a_s_p*p.Area*p.L_p);

%%% SOLID PHASE DYNAMICS
% Solid phase diffusivity temperature dependence
D_s_n0 = anodeD(SOC_n, p);
D_s_p0 = cathodeD(SOC_p, p);

p.D_s_n = D_s_n0; 
p.D_s_p = D_s_p0;

% Construct (A,B) matrices for solid-phase Li diffusion
[A_n,A_p,B_n,B_p,C_n,C_p,D_n,D_p] = spm_plant_obs_mats(p);

% Compute surface concentrations
c_ss_n = C_n*c_s_n + D_n*jn_tot;
c_ss_p = C_p*c_s_p + D_p*jp_tot;
% Remark: I am cheating slightly here. jn_tot should be jn, but doing so
% imposes an algebraic equation. This forms a DAE. I am going to
% approximate jn by jn_tot, which should be ok, since jn and jn_tot have a
% error on the order of 0.001%

%%% ELECTROLYTE PHASE DYNAMICS
% Electrolyte error handling
% if any(c_e<0)
%     c_e0 = evalin('base','ce0');
%     c_e(1:(p.Nxn-1)) = c_e0(1:(p.Nxn-1));
% end

% Compute electrolyte Boundary Conditions
c_e_bcs = p.ce.C * c_e;

ce0n = c_e_bcs(1);
cens = c_e_bcs(2);
cesp = c_e_bcs(3);
ce0p = c_e_bcs(4);

% Electrolyte error handling
% elektrolyte_flag = 0;
% if (c_e<0 | ce0n<0 | ce0p<0) 
%     elektrolyte_flag = 1;
%     ce0n = 0.5;
% end

% Separate and aggregate electrolyte concentration
c_en = c_e(1:(p.Nxn-1));
c_es = c_e((p.Nxn-1)+1:(p.Nxn-1)+(p.Nxs-1));
c_ep = c_e((p.Nxn-1)+p.Nxs : end);
c_ex = [ce0n; c_en; cens; c_es; cesp; c_ep; ce0p];


%% Voltage output

% Average electrolyte concentrations
cen_bar = mean(c_ex(1:p.Nxn+1,:));
ces_bar = mean(c_ex((p.Nxn+1):(p.Nxn+p.Nxs+1),:));
cep_bar = mean(c_ex((p.Nxn+p.Nxs+1):(p.Nxn+p.Nxs+p.Nxp+1),:));

% Overpotentials due to electrolyte subsystem
kap_n = electrolyteCond(cen_bar,T1,p);
kap_s = electrolyteCond(ces_bar,T1,p);
kap_p = electrolyteCond(cep_bar,T1,p);

% Bruggeman relationships
kap_n_eff = kap_n * p.epsilon_e_n.^(p.brug);
kap_s_eff = kap_s * p.epsilon_e_s.^(p.brug);
kap_p_eff = kap_p * p.epsilon_e_p.^(p.brug);

% Activity coefficient
dfca_n = electrolyteAct(cen_bar);
dfca_s = electrolyteAct(ces_bar);
dfca_p = electrolyteAct(cep_bar);

% Kinetic reaction rate, adjusted for Arrhenius temperature dependence
k_n0 = anodeK(SOC_n, p);
k_p0 = cathodeK(SOC_p, p);

p.k_n = k_n0;
p.k_p = k_p0;

% Stochiometric Concentration Ratio
theta_n = c_ss_n / p.c_s_n_max;
theta_p = c_ss_p / p.c_s_p_max;

% Equilibrium Potential
Unref = refPotentialAnode(p,theta_n);
Upref = refPotentialCathode(p,theta_p);

% Exchange current density
c_e_bar = [cen_bar; ces_bar; cep_bar];
[i_0n,i_0p] = exch_cur_dens(p,c_ss_n,c_ss_p,c_e_bar);

% Overpotentials
RTaF=(p.R*T1)/(p.alph*p.Faraday);
eta_n = RTaF * asinh(-cur / (2*p.a_s_n*p.Area*p.L_n*i_0n(1)));
eta_p = RTaF * asinh(cur / (2*p.a_s_p*p.Area*p.L_p*i_0p(end)));

% Total resistance (film + growing SEI layer)
R_tot_n = p.R_f_n + delta_sei/p.kappa_P;
R_tot_p = p.R_f_p + 0;

% SPM Voltage (i.e. w/o electrolyte concentration terms)
eta_sei_n = (R_tot_n/(p.a_s_n*p.L_n*p.Area))*-cur;
eta_sei_p = (R_tot_p/(p.a_s_p*p.L_p*p.Area))*-cur;

V_noVCE = eta_p - eta_n + Upref - Unref - eta_sei_n + eta_sei_p;

% Overpotential due to electrolyte conductivity
V_electrolyteCond = (p.L_n/(2*kap_n_eff) + 2*p.L_s/(2*kap_s_eff) + p.L_p/(2*kap_p_eff))*(cur/p.Area);

% Overpotential due to electrolyte polarization
V_electrolytePolar = (2*p.R*T1)/(p.Faraday) * (1-p.t_plus)* ...
        ( (1+dfca_n) * (log(cens) - log(ce0n)) ...
         +(1+dfca_s) * (log(cesp) - log(cens)) ...
         +(1+dfca_p) * (log(ce0p) - log(cesp)));

% Add 'em up!
V = V_noVCE + V_electrolyteCond + V_electrolytePolar;


%% Aging Dynamics

% Molar flux of side rxn [mol/s-m^2]
j_s = 0;

% Molar flux of intercalation
jn = (abs(jn_tot) - abs(j_s)) * sign(jn_tot);
jp = jp_tot;


%% Solid Phase Dynamics

% ODE for c_s
c_s_n_dot = A_n*c_s_n + B_n*jn;
c_s_p_dot = A_p*c_s_p + B_p*jp;


%% Electrolyte Dynamics

% Compute Electrolyte Diffusion Coefficient and Derivative
[D_en0,dD_en0] = electrolyteDe(c_en);
[D_es0,dD_es0] = electrolyteDe(c_es);
[D_ep0,dD_ep0] = electrolyteDe(c_ep);

% Adjustment for Arrhenius temperature dependence
Arrh_De = exp(p.E.De/p.R*(1/p.T_ref_electrolyte - 1/T1));
D_en = D_en0 * Arrh_De;
D_es = D_es0 * Arrh_De;
D_ep = D_ep0 * Arrh_De;
dD_en = dD_en0 * Arrh_De;
dD_es = dD_es0 * Arrh_De;
dD_ep = dD_ep0 * Arrh_De;

% Apply bruggeman relation 
D_en_eff = D_en .* p.epsilon_e_n.^(p.brug-1);
dD_en_eff = dD_en .* p.epsilon_e_n.^(p.brug-1);

D_es_eff = D_es .* p.epsilon_e_s.^(p.brug-1);
dD_es_eff = dD_es .* p.epsilon_e_s.^(p.brug-1);

D_ep_eff = D_ep .* p.epsilon_e_p.^(p.brug-1);
dD_ep_eff = dD_ep .* p.epsilon_e_p.^(p.brug-1);

% System Matrices have all been precomputed & stored in param struct "p"

% Compute derivative
c_en_dot = dD_en_eff.*(p.ce.M1n*c_en + p.ce.M2n*c_e_bcs(1:2)).^2 ...
    + D_en_eff.*(p.ce.M3n*c_en + p.ce.M4n*c_e_bcs(1:2)) + diag(p.ce.M5n)*jn;

c_es_dot = dD_es_eff.*(p.ce.M1s*c_es + p.ce.M2s*c_e_bcs(2:3)).^2 ...
    + D_es_eff.*(p.ce.M3s*c_es + p.ce.M4s*c_e_bcs(2:3));

c_ep_dot = dD_ep_eff.*(p.ce.M1p*c_ep + p.ce.M2p*c_e_bcs(3:4)).^2 ...
    + D_ep_eff.*(p.ce.M3p*c_ep + p.ce.M4p*c_e_bcs(3:4)) + diag(p.ce.M5p)*jp;

% Assemble c_e_dot
c_e_dot = [c_en_dot; c_es_dot; c_ep_dot];

% Electrolyte error handling
% if  elektrolyte_flag==1 & c_en_dot(1)<0
%     c_e_dot = zeros(size(c_e_dot));
% end

%% Thermal Dynamics

% State-of-Charge (Bulk)
r_vec = (0:p.delta_r_n:1)';
c_n = [c_s_n(1); c_s_n; c_ss_n];
c_p = [c_s_p(1); c_s_p; c_ss_p];
SOC_n = 3/p.c_s_n_max * trapz(r_vec,r_vec.^2.*c_n);
SOC_p = 3/p.c_s_p_max * trapz(r_vec,r_vec.^2.*c_p);

%% Concatenate time derivatives
x_dot = [c_s_n_dot; c_s_p_dot; c_e_dot];

%% Electrolyte error handling
% assignin('base','ce0',c_e);

%% Concatenate outputs
varargout{1} = V;
varargout{2} = V_noVCE;
varargout{3} = SOC_n;
varargout{4} = SOC_p;
varargout{5} = c_ss_n;
varargout{6} = c_ss_p;
varargout{7} = c_ex';
varargout{8} = eta_n;
varargout{9} = eta_p;
varargout{10} = Unref;
varargout{11} = Upref;
varargout{12} = V_electrolytePolar;
varargout{13} = V_electrolyteCond;  
varargout{14} = eta_sei_n; 
varargout{15} = eta_sei_p; 

