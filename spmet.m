%% Single Particle Model w/ Electrolyte & Temperature
%   Originally published by Scott Moura under GPL-3.0 license in December, 2016.
%   Link to original repository: https://github.com/scott-moura/SPMeT
%   Restructured and extended by Nikolaos Wassiliadis, Andreas Bach, Matthias Wanzel in April, 2020.

addpath(genpath(pwd));
clear;
clc;
close all;

%% Electrochemical model parameter
% Load US18650VTC5A NCA cell parameter
run parameter/parameter.m

% Load measurement
load('US18650VTC5A_1C_20deg.mat');

% Calculate electrode capacities
[cn_low,cp_low] = init_SOC(p,0); % Calculate concentrations at 0% SOC
[cn_high,cp_high] = init_SOC(p,100);  % Calculate concentrations at 100% SOC
Delta_cn = cn_high-cn_low; % Calculate anode concentration difference
Delta_cp = cp_low-cp_high; % Calculate cathode concentration difference
qn = p.epsilon_s_n*p.Area*p.L_n*Delta_cn*p.Faraday/3600; % Calculate anode capacity
qp = p.epsilon_s_p*p.Area*p.L_p*Delta_cp*p.Faraday/3600; % Calculate cathode capacity

% Cell capacity deviation
Dataset.Q = 2.52; % Actual full cell capacity
Q_t = Dataset.Q/qn; % Calculate deviation to original parameterization
p.Q_n = p.Q_n*Q_t; % Scale anode capacity accordingly
p.Q_p = p.Q_p*Q_t;  % Scale cathode capacity accordingly
p.L_n = p.Q_n*3600/(p.epsilon_s_n*p.Area*p.c_s_n_max*p.Faraday); % Adjust anode thickness
p.L_p = p.Q_p*3600/(p.epsilon_s_p*p.Area*p.c_s_p_max*p.Faraday); % Adjust cathode thickness

% Interpolate data to seconds
Dataset.Time = Dataset.Time*3600; % Convert time to seconds
ip_t = 0:1:max(Dataset.Time); % Normalize time to 1s distance
Dataset.U = interp1(Dataset.Time, Dataset.U, ip_t); % Interpolate voltage
Dataset.I = interp1(Dataset.Time, Dataset.I, ip_t); % Interpolate current
Dataset.T2 = interp1(Dataset.Time, Dataset.T1, ip_t); % Interpolate measured temperature

% Parse data
t = ip_t; % Assign time as input
I =  Dataset.I; % Assign current as input
T2 = Dataset.T2; % Assing temperature as input

%% Input Data structure with time, current and temperature
% Current | Positive = Charge | Negative = Discharge
data.time = t;          % [s]
data.cur = I;           % [A]
data.T2 = T2+273.15;    % [K]

%% Preallocation and initial conditions

% Finite difference for spherical particle
p.Nr = 20; % Particle resolution
Nr = p.Nr; % Assign particle resolution as input
p.delta_r_n = 1/p.Nr; % Anode fraction percentage 
p.delta_r_p = 1/p.Nr; % Cathode fraction percentage
r_vec = (0:p.delta_r_n:1)'; % Electrode fraction vector
r_vecx = r_vec(2:end-1); % Realign vector

% Finite difference points along x-coordinate
p.Nxn = 10; % Anode resolution
p.Nxs = 5; % Seperator resolution
p.Nxp = 10; % Cathode resolution
p.Nx = p.Nxn+p.Nxs+p.Nxp; % Axial resolution
Nx = p.Nx - 3; % Assign axial resolution as input
x_vec_spme = linspace(0,1,Nx+4); % Vector for solution space

p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;


%% Set initial conditions

% Solid concentration
SOC_0 = 0; % Set inital full cell SOC
[csn0,csp0] = init_SOC(p,SOC_0); % Calculate initial solid concentration
c_n0 = csn0 * ones(p.Nr-1,1); % Preallocate
c_p0 = csp0 * ones(p.Nr-1,1); % Preallocate

% Electrolyte concentration
ce0 = p.c_e*ones(Nx,1); % Preallocate

% Temperature
T10 = data.T2(1); % Initial cell temperature

% SEI layer
delta_sei0 = 0; % Initial SEI layer


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
disp('Simulating SPMeT...');

% Initial Conditions
x0 = [c_n0; c_p0; ce0; T10; delta_sei0]; % Set initial conditions for simulation
options = odeset('Stats','off','MaxStep',10); % Set solver settings

% Integrate
[t,x] = ode23s(@(t,x) ode_spmet(t,x,data,p),t,x0,options); % Solve SPMeT

% Parse states
c_s_n = x(:,1:(p.Nr-1));
c_s_p = x(:,p.Nr : 2*(p.Nr-1));
c_ex = x(:,2*p.Nr-1 : 2*p.Nr-1+p.Nx-4);
T1 = x(:,end-1);
delta_sei = x(:,end);

%% Output function
NT = length(data.time); % Assign all sample points to counter variable NT

V = zeros(NT,1);                        % Voltage SPMeT
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
    [~,V(k),V_spm(k),SOC_n(k),SOC_p(k),c_ss_n(k),c_ss_p(k),c_e(:,k),eta_n(k),eta_p(k),Unref(k),Upref(k),V_ep(k),...
        V_ec(k), eta_sei_n(k), eta_sei_p(k)] = ode_spmet(t(k),x(k,:)',data,p);

    % Aggregate solid concentrations
    c_n(k,:) = [c_s_n(k,1), c_s_n(k,:), c_ss_n(k)];
    c_p(k,:) = [c_s_p(k,1), c_s_p(k,:), c_ss_p(k)];
    
    % Total moles of lithium in solid
    n_Li_s(k) = (3*p.epsilon_s_p*p.L_p*p.Area) * trapz(r_vec,r_vec.^2.*c_p(k,:)') ...
            + (3*p.epsilon_s_n*p.L_n*p.Area) * trapz(r_vec,r_vec.^2.*c_n(k,:)');
        
    if k > 1
        fprintf(reverseStr);
    end
    msg = sprintf('Processed %d/%d', k, NT);
    fprintf(msg);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

% Scale to system
V = V.*p.bat_s + p.R_ECR/p.bat_p *2 *I' *p.bat_s;

% Output elapsed simulation time
simtime = toc;
fprintf(1,'\n Elapsed time: %4.1f sec or %2.2f min \n',simtime,simtime/60);

% Visualize results
figure;
plot(t, V, 'b')
hold on;
plot(t, Dataset.U, 'k')
hold on;
legend({'$$U_{Simulation}$$','$$U_{Measurement}$$'},'interpreter','latex','Fontsize',12,'Location','best');
ylabel('U in V');
xlabel('Time in s');
