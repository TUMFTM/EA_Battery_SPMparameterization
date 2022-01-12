%% Single Particle Model w/ Electrolyte & Temperature as a function
%   See spmet.m for credentials

function [U_sim] = spmet_fnc(factors, t, I, T, SOC0, Nr, Nxn, Nxs, Nxp, max_step)
% factors = factors to scale regression functions of k [anode cathode]
% t = Time vector in h
% I = Current vector in A
% T = Ambient temperature scalar in K (iso-thermal)
% CA = Actual cell capacity in Ah
% Nr = Radial FDM discretization
% Nxn = Axial FDM dicreztization in anode
% Nxs = Axial FDM dicreztization in separator
% Nxp = Axial FDM dicreztization in cathode

% Allow only positive factors for full cell level scaling
if isempty(find(factors < 0, 4))
    %% Electrochemical model parameter
    % Load US18650VTC5A NCA cell parameter
    run parameter/parameter.m

    % Set ambient temperature to default
    p.T_amb = T;

    % Set factors to scale regression function of activation rate
    p.K_multi_n = factors(1,1);
    p.K_multi_p = factors(1,2);

    % Set factors to scale regression function of diffusion coefficient
    p.D_multi_n = factors(1,3);
    p.D_multi_p = factors(1,4);

    %% Input Data structure with time, current and temperature
    % Current | Positive <=> Charge | Negative <=> Discharge
    data.time = t;      % [s]
    data.cur = I;       % [A]

    %% Preallocation and initial conditions

    % Finite difference for spherical particle
    p.Nr = Nr; % Particle resolution
    % Nr = p.Nr; % Hand-over
    p.delta_r_n = 1/p.Nr; % Anode fraction percentage 
    p.delta_r_p = 1/p.Nr; % Cathode fraction percentage
    % r_vec = (0:p.delta_r_n:1)'; % Electrode fraction vector
    % r_vecx = r_vec(2:end-1);

    % Finite difference points along x-coordinate
    p.Nxn = Nxn; % Anode resolution
    p.Nxs = Nxs; % Seperator resolution
    p.Nxp = Nxp; % Cathode resolution
    p.Nx = p.Nxn+p.Nxs+p.Nxp; % Axial resolution
    Nx = p.Nx - 3;
    % x_vec_spme = linspace(0,1,Nx+4);

    p.delta_x_n = 1 / p.Nxn;
    p.delta_x_s = 1 / p.Nxs;
    p.delta_x_p = 1 / p.Nxp;

    %% Set initial conditions

    % Solid concentration
    SOC_0 = SOC0*100;
    [csn0,csp0] = init_SOC(p,SOC_0); % Calculate initial solid concentration
    c_n0 = csn0 * ones(p.Nr-1,1); % Preallocate
    c_p0 = csp0 * ones(p.Nr-1,1); % Preallocate

    % Electrolyte concentration
    ce0 = p.c_e*ones(Nx,1); % Preallocate

    % Temperature
    % T10 = data.T2(1); % Initial cell temperature

    % SEI layer
    % delta_sei0 = 0; % Initial SEI layer

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
    tic; % Start time tracking

    % Initial Conditions
    x0 = [c_n0; c_p0; ce0]; % Set initial conditions for simulation
    options = odeset('MaxStep', max_step,...
                     'Stats','off'); % See solver results

    % Integrate
    [t_sim,x] = ode23s(@(t,x) ode_spmet(t,x,data,p),t,x0,options); % Solve SPMeT

    % Parse states
    c_s_n = x(:,1:(p.Nr-1));
    c_s_p = x(:,p.Nr : 2*(p.Nr-1));
    % c_ex = x(:,2*p.Nr-1 : 2*p.Nr-1+p.Nx-4);
    % T1 = x(:,end-1);
    % delta_sei = x(:,end);

    %% Output function
    NT = length(data.time); % Assign all sample points to counter variable NT

    U_sim = zeros(NT,1);                    % Voltage SPMeT
    V_spm = zeros(NT,1);                    % Voltage without electrolyte 
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
        [~,U_sim(k),V_spm(k),SOC_n(k),SOC_p(k),c_ss_n(k),c_ss_p(k),c_e(:,k),eta_n(k),eta_p(k),Unref(k),Upref(k),V_ep(k),...
            V_ec(k), eta_sei_n(k), eta_sei_p(k)] = ode_spmet(t_sim(k),x(k,:)',data,p); 
    end

    % Output elapsed simulation time
    simtime = toc;
    
else
    % if factors are negative, set simulation voltage to zero
    NT    = length(t);
    U_sim = zeros(NT,1);
end

U_sim = real(U_sim);

end
