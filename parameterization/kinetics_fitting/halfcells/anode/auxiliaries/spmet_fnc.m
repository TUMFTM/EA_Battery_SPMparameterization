%% Single Particle Model w/ Electrolyte & Temperature as a function
%   See spmet.m for credentials

function [t_sim, U_sim, varargout] = spmet_fnc(t, I, p0, SOC0, Nr, Nxn, Nxs, Nxp)
% Inputs:
% t = Time vector in s
% I = Current vector in A
% p0 = Default parameter for electrode diffusion and reaction rate in [m^2/s], [(A/m^2)*(mol^3/mol)^(1+alpha)]
% SOC0 = Initial State-of-Charge in -
% Nr = Radial FDM discretization in -
% Nxn = Axial FDM dicreztization in anode in -
% Nxs = Axial FDM dicreztization in separator in -
% Nxp = Axial FDM dicreztization in cathode in -
% Outputs:
% t_sim = Simulation time vector in s
% U_sim = Simulation voltage vector in V
% varargout = Variable outputs, see end of function

% Allow only positive factors for full cell level scaling
if isempty(find(p0 < 0, 2)) %&& isempty(find(p0 > 5, 2))

    %% Electrochemical model parameter
    % Load US18650VTC5A NCA cell parameter
    run parameter/parameter.m

    % Extend parameter with default concentrations and reaction rates
    p.D_s_n0 = 7e-15*p0(1);
    p.D_s_p0 = 7e-15;
    p.k_n0 = 3e-6*p0(2);
    p.k_p0 = 3e-6;

    %% Input Data structure with time, current and temperature
    % Current | Positive <=> Charge | Negative <=> Discharge
    data.time = t;      % [s]
    data.cur = I;       % [A]

    %% Preallocation and initial conditions

    % Finite difference for spherical particle
    p.Nr = Nr; % Particle resolution
    %Nr = p.Nr; % Hand-over
    p.delta_r_n = 1/p.Nr; % Anode fraction percentage 
    p.delta_r_p = 1/p.Nr; % Cathode fraction percentage
    %r_vec = (0:p.delta_r_n:1)'; % Electrode fraction vector
    %r_vecx = r_vec(2:end-1);

    % Finite difference points along x-coordinate
    p.Nxn = Nxn; % Anode resolution
    p.Nxs = Nxs; % Seperator resolution
    p.Nxp = Nxp; % Cathode resolution
    p.Nx = p.Nxn+p.Nxs+p.Nxp; % Axial resolution
    Nx = p.Nx - 3;
    %x_vec_spme = linspace(0,1,Nx+4);

    p.delta_x_n = 1 / p.Nxn;
    p.delta_x_s = 1 / p.Nxs;
    p.delta_x_p = 1 / p.Nxp;

    %% Set initial conditions

    % Solid concentration
    SOC0 = SOC0*100;
    [csn0,csp0] = init_SOC(p,SOC0);
    c_n0 = csn0 * ones(p.Nr-1,1); % Preallocate
    c_p0 = csp0 * ones(p.Nr-1,1); % Preallocate

    % Electrolyte concentration
    ce0 = p.c_e*ones(Nx,1); % Preallocate

    % Temperature
    %T10 = p.T_amb; % Iso-thermal setup

    % SEI layer
    %delta_sei0 = 0; % Initial SEI layer

    %% Generate constant system matrices

    % Electrolyte concentration matrices
    [M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p);

    p.ce.M1n = M1n;
    p.ce.M2n = M2n;
    p.ce.M3n = M3n;
    p.ce.M4n = M4n;
    p.ce.M5n = M5n;

    p.ce.M1s = M1s;
    p.ce.M2s = M2s;
    p.ce.M3s = M3s;
    p.ce.M4s = M4s;

    p.ce.M1p = M1p;
    p.ce.M2p = M2p;
    p.ce.M3p = M3p;
    p.ce.M4p = M4p;
    p.ce.M5p = M5p;

    p.ce.C = C_ce;

    clear M1n M2n M3n M4n M5n M1s M2s M3s M4s M1p M2p M3p M4p M5p C_ce;

    %% Simulate SPMeT plant

    % Initial Conditions
    x0 = [c_n0; c_p0; ce0]; % Set initial conditions for simulation
    options = odeset('MaxStep',100, ... % Relax error tolerances to allow for easier ODE solving
                     'Stats','off'); 	% See solver results

    % Integrate
    [t,x] = ode23s(@(t,x) ode_spmet(t,x,data,p),t,x0,options); % Solve SPMeT

    % Parse states
    c_s_n = x(:,1:(p.Nr-1));
    c_s_p = x(:,p.Nr : 2*(p.Nr-1));
    %c_ex = x(:,2*p.Nr-1 : 2*p.Nr-1+p.Nx-4);
    %T1 = x(:,end-1);
    %delta_sei = x(:,end);

    %% Output function
    NT = length(data.time); % Assign all sample points to counter variable NT

    t_sim = t;                              % Simulation time
    U_sim = zeros(NT,1);                    % Voltage SPMeT
    U_spm = zeros(NT,1);                    % Voltage without electrolyte 
    SOC_n = zeros(NT,1);                    % State of charge neg. electrode 
    SOC_p = zeros(NT,1);                    % State of charge pos. electrode 
    c_ss_n = zeros(NT,1);                   % Surface concentration neg. particle
    c_ss_p = zeros(NT,1);                   % Surface concentration pos. particle
    c_n = zeros(NT,p.Nr+1);                 % Solid concentration neg. electrode
    c_p = zeros(NT,p.Nr+1);                 % Solid concentration pos. electrode
    c_e = zeros(p.Nx+1,NT);                 % Electrolyte concentrations
    n_Li_s = zeros(NT,1);                   % Total Moles of Lithium in Solid
    eta_n = zeros(NT,1);                    % Activation overpotential neg. electrode based on Butler–Volmer equation
    eta_p = zeros(NT,1);                    % Activation overpotential pos. electrode based on Butler–Volmer equation
    Unref = zeros(NT,1);                    % Equilibrium potential neg. electrode
    Upref = zeros(NT,1);                    % Equilibrium potential pos. electrode
    V_ep = zeros(NT,1);                     % Overpotential due to electrolyte polarization
    V_ec = zeros(NT,1);                     % Overpotential due to electrolyte conductivity
    eta_sei_n = zeros(NT,1);                % Overpotential neg. electrode due to SEI resistance
    eta_sei_p = zeros(NT,1);                % Overpotential pos. electrode due to SEI resistance

    for k = 1:NT
        % Compute outputs
        [~,U_sim(k),U_spm(k),SOC_n(k),SOC_p(k),c_ss_n(k),c_ss_p(k),c_e(:,k),eta_n(k),eta_p(k),Unref(k),Upref(k),V_ep(k),...
            V_ec(k), eta_sei_n(k), eta_sei_p(k)] = ode_spmet(t(k),x(k,:)',data,p);        
    end

    % Assign auxiliaries to output
    varargout{1} = U_spm;
    varargout{2} = SOC_n;
    varargout{3} = SOC_p;
    varargout{4} = c_ss_n;
    varargout{5} = c_ss_p;
    varargout{6} = c_e;
    varargout{7} = eta_n;
    varargout{8} = eta_p;
    varargout{9} = Unref;
    varargout{10} = Upref;
    varargout{11} = eta_sei_n;
    varargout{12} = eta_sei_p;
    varargout{13} = c_s_n;
    varargout{14} = c_s_p;

else % if factors are negative, set simulation voltage to zero
    %% Output function
    NT = length(t); % Assign all sample points to counter variable NT
    
    t_sim = zeros(NT,1);                    % Simulation time
    U_sim = zeros(NT,1);                    % Voltage SPMeT
    
    % Assign auxiliaries to output
    varargout{1} = 0;
    varargout{2} = 0;
    varargout{3} = 0;
    varargout{4} = 0;
    varargout{5} = 0;
    varargout{6} = 0;
    varargout{7} = 0;
    varargout{8} = 0;
    varargout{9} = 0;
    varargout{10} = 0;
    varargout{11} = 0;
    varargout{12} = 0;
    varargout{13} = 0;
    varargout{14} = 0;
end
end
