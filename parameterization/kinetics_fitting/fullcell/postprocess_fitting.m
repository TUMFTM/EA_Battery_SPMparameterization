%% Postprocessing of optimized factors used for Arrhenius fitting
%  Developed by Nikolaos Wassiliadis and Andreas Bach in November, 2020

% addpath(genpath(pwd));
addpath('./auxiliaries');
addpath('./measurement');
addpath('./parameter');
clear;
close all;
clc;

%% Get measurement data

temp263 = load('intermediate_results_factors_263_Kelvin');
temp278 = load('intermediate_results_factors_278_Kelvin');
temp293 = load('intermediate_results_factors_293_Kelvin');
temp308 = load('intermediate_results_factors_308_Kelvin');
temp323 = load('intermediate_results_factors_323_Kelvin');

%% Prepare and execute Arrhenius k fitting

% x values used in k & D fitting including anode and cathode
x     = [-10, 5, 20, 35, 50] + 273.15;

% y values anode k
y_k_n = [temp263.x_opt(1), temp278.x_opt(1), temp293.x_opt(1), temp308.x_opt(1), temp323.x_opt(1)];
% Arrhenius fitting anode
[coeff_k_n] = arrh_createFit(x, y_k_n, 'k_n');

% y values cathode k
y_k_p = [temp263.x_opt(2), temp278.x_opt(2), temp293.x_opt(2), temp308.x_opt(2), temp323.x_opt(2)];
% Arrhenius fitting cathode
[coeff_k_p] = arrh_createFit(x, y_k_p, 'k_p');

%% Prepare and execute Arrhenius D fitting

% y values anode D
y_D_n = [temp263.x_opt(3), temp278.x_opt(3), temp293.x_opt(3), temp308.x_opt(3), temp323.x_opt(3)];
% Arrhenius fitting anode
[coeff_D_n] = arrh_createFit(x, y_D_n, 'D_n');

% y values cathode D
y_D_p = [temp263.x_opt(4), temp278.x_opt(4), temp293.x_opt(4), temp308.x_opt(4), temp323.x_opt(4)];
% Arrhenius fitting cathode
[coeff_D_p] = arrh_createFit(x, y_D_p, 'D_p');

%% Save results

% Isolate factors at reference temperature 293.15K
ref_factors         = [coeff_k_n(1), coeff_k_p(1), coeff_D_n(1), coeff_D_p(1)];

% Isolate activation energies at reference temperature 293.15K
activation_energies = [coeff_k_n(2), coeff_k_p(2), coeff_D_n(2), coeff_D_p(2)];

save('kinetics_fullcell', 'ref_factors', 'activation_energies')
savefig('final_results_arrhenius')