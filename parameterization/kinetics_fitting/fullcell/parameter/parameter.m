%% Electrochemical cell parameter
%   Originally published by Scott Moura under GPL-3.0 license in December, 2016.
%   Link to original repository: https://github.com/scott-moura/SPMeT
%   Reorganized and extended with a parameter set of a commercially available US18650VTC5A NCA
%   lithium-ion cell by Nikolaos Wassiliadis, Andreas Bach, Matthias Wanzel in April, 2020.

% Choose cell capacity depending on the measurement temperature
switch T
    case 263.15
        Q_cell = 2.324876864349120;
    case 278.15
        Q_cell = 2.383004312947400;
    case 293.15
        Q_cell = 2.417569095185370;
    case 308.15
        Q_cell = 2.497250669594920;
    case 323.15
        Q_cell = 2.504991235060170;
end

%% Load OCP fitting parameter

% Balancing and Alignment 
load('parameter/B_A.mat');
p.B_A = Results;
% Anode OCV 
transfer = load('parameter/OCV_anode.mat');
p.p_n = transfer.p;
% Cathode OCV 
transfer = load('parameter/OCV_cathode.mat');
p.p_p = transfer.p;
% Anode Kinetics
transfer = load('parameter/kinetics_anode');
p.kin_k_n = transfer.anode_fitk;
p.kin_D_n = transfer.anode_fitD;
% Cathode Kinetics
transfer = load('parameter/kinetics_cathode');
p.kin_k_p = transfer.cathode_fitk;
p.kin_D_p = transfer.cathode_fitD;

%% Load geometric parameter

p.Area = 0.1024;  % Electrode area [m^2] - Lain (2019) 

% Thickness of layers
p.L_s = 8e-6;     % Lain (2019) Thickness of separator [m]
L_ccn = 14e-6;    % Lain (2019) Thickness of negative current collector [m]
L_ccp = 15e-6;    % Lain (2019) Thickness of positive current collector [m]

% Particle Radii
p.R_s_n = 7e-6;   % Lain (2019) Radius of solid particles in negative electrode [m]
p.R_s_p = 2.5e-6; % Lain (2019) Radius of solid particles in positive electrode [m]

% Volume fractions
p.epsilon_s_n = 0.73*0.95;      % Volume fraction in solid for neg. electrode - Lain (2019) 95% active material 
p.epsilon_s_p = 0.87*0.96;      % Volume fraction in solid for pos. electrode - Lain (2019) 96% active material

p.epsilon_e_n = 0.27;   % Lain (2019) Volume fraction in electrolyte for neg. electrode
p.epsilon_e_s = 0.4;    % Volume fraction in electrolyte for separator
p.epsilon_e_p = 0.13;   % Lain (2019) Volume fraction in electrolyte for pos. electrode

% make element to caclulate phi_{s}
p.epsilon_f_n = 1-p.epsilon_s_n-p.epsilon_e_n;  % Volume fraction of filler in neg. electrode
p.epsilon_f_p = 1-p.epsilon_s_p-p.epsilon_e_p;  % Volume fraction of filler in pos. electrode

epsilon_f_n = p.epsilon_f_n;  % Volume fraction of filler in neg. electrode
epsilon_f_p = p.epsilon_f_p;  % Volume fraction of filler in pos. electrode

% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]

% Mass densities
rho_sn = 2250;    % CAS 7782-42-5 2200 - Lain (2019) 2250
rho_sp = 4800;    % CAS 193214-24-3 4450 - Lain (2019) 4800
rho_e =  1324;    % Electrolyte [kg/m^3]
rho_f_n = 1743;   % Filler negativ [kg/m^3] Lain (2019) Binder:Carbon 4:1  
rho_f_p = 1812;   % Filler positive [kg/m^3] Lain (2019) Binder:Carbon 1:1 
rho_ccn = 8954;   % Current collector in negative electrode
rho_ccp = 2707;   % Current collector in positive electrode


%% Transport parameter

p.brug = 1.5;       % Bruggeman porosity

% Miscellaneous
p.t_plus = 0.38;            % Transference number - Valoen (2005) 
p.Faraday = 96485.33289;    % Faraday's constant, [Coulumbs/mol]


%% Kinetic parameter

p.R = 8.314472;         % Gas constant, [J/mol-K]

p.alph = 0.5;           % Charge transfer coefficients

p.R_f_n = 4.4e-3;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_p = 0;            % Resistivity of SEI layer, [Ohms*m^2]


%% Aging submodel Params (Unused)

%   SEI Layer Growth model
%   Adopted from Ramadass et al (2004) [Univ of South Carolina]
%   "Development of First Principles Capacity Fade Model for Li-Ion Cells"
%   DOI: 10.1149/1.1634273
%   NOTE: These parameters have NOT been experimentally validated by eCAL

p.kappa_P = 1;      % [S/m] conductivity of side rxn product
p.M_P = 7.3e1;      % [kg/mol] molecular weight of side rxn product
p.rho_P = 2.1e3;    % [kg/m^3] mass density of side rxn product
p.i0s = 0; %1.5e-6; % [A/m^2] exchange current density of side rxn
p.Us = 0.4;         % [V] reference potential of side rxn


%% Concentrations

q_n = 343;                                         % Active material capacity, [mAh/g] - Lain (2019)
p.c_s_n_max = 3.6e3 * q_n * rho_sn / p.Faraday;    % Max concentration in anode, [mol/m^3]

q_p = 168;                                         % Active material capacity, [mAh/g] - Lain (2019)
p.c_s_p_max = 3.6e3 * q_p * rho_sp / p.Faraday;    % Max concentration in cathode, [mol/m^3]

p.c_e = 1e3;                                       % Fixed electrolyte concentration, [mol/m^3]


%% Cell Deviation
Q_t = Q_cell/2.581737877888906; % Capacity tuning Parameter  

% Match thickness of electrodes according to results from Balancing and Alignment
p.Q_n = (max(p.B_A.Q_anode) - min(p.B_A.Q_anode))/3600*Q_t;      % Anode capacity, [Ah]
p.Q_p = (max(p.B_A.Q_cathode) - min(p.B_A.Q_cathode))/3600*Q_t;  % Cathode capacity according, [Ah]
p.L_n = p.Q_n*3600/(p.epsilon_s_n*p.Area*p.c_s_n_max*p.Faraday); % Thickness of neg. electrode, [m]
p.L_p = p.Q_p*3600/(p.epsilon_s_p*p.Area*p.c_s_p_max*p.Faraday); % Thickness of pos. electrode, [m]


%% Lumped weight of cell 

% Compute cell mass [kg/m^2]
m_n = p.L_n * (rho_e*p.epsilon_e_n + rho_sn*p.epsilon_s_n + rho_f_n*epsilon_f_n);
m_s = p.L_s * (rho_e*p.epsilon_e_s);
m_p = p.L_p * (rho_e*p.epsilon_e_p + rho_sp*p.epsilon_s_p + rho_f_p*epsilon_f_p);
m_cc = rho_ccn*L_ccn + rho_ccp*L_ccp;

% Lumped density [kg/m^2]
p.rho_avg = m_n + m_s + m_p + m_cc;

% Lumped weight [kg]
p.w_cell = p.rho_avg*p.Area;


%% Thermodynamic parameter

% Reference temperature
p.T_ref = 293.15;               % Arrhenius reference temperature [K]
p.T_ref_electrolyte = 296.15;   % Arrhenius reference temperature for electrolyte [K]

% Activation Energies
p.E.De = 17.12e3;               % Diffusion coeff electrolyte [J/mol]
p.E.kappa_e =17.12e3;           % Activation energy electrolyte conductivity [J/mol]

% Heat transfer parameters
p.C1 = 986.2*p.w_cell;          % Heat capacity [J/K]
p.h12 = 0.0306;                 % Heat transfer coefficient [W/K] alpha = 3.4 [W/(m^2*K)]