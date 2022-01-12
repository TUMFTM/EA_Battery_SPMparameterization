%% Create regression of half cell current step response fitting
%  Developed by Nikolaos Wassiliadis, Andreas Bach and Matthias Wanzel in November, 2020

addpath('./measurement');
addpath('./parameter');
addpath('./auxiliaries');
clear;
clc;
close all;

%% Load fitting (reaction rate and diffusion coefficient)
load('intermediate_results_coefficients_v5');
GITT = load('Halfcell_GITT');

%% Fit regression

% Prepare curve fitting
[~,D] = prepareCurveData(coeff(:,1),coeff(:,2)); % Diffusion coefficients
[x,k] = prepareCurveData(coeff(:,1),coeff(:,3)); % Reaction rate

% Exclude errorenous data points, if necessary
exD = excludedata(x, D, 'Indices', find(D<7e-16));
exk = excludedata(x, k, 'Indices', [2 4:5 8:10 38]);

% Define fitting functions
fD = fittype('poly3');
fk = fittype('poly2'); % Reaction rate

% Fit
opts = fitoptions('Method', 'NonlinearLeastSquares', 'Algorithm', 'levenberg-marquardt');
opts.Exclude = exD; % Excluded points if necessary
[anode_fitD, ~] = fit(x, D, fD, opts);
opts = fitoptions( 'Method', 'NonlinearLeastSquares','Algorithm', 'levenberg-marquardt');
opts.Exclude = exk; % Excluded points if necessary
[anode_fitk, ~] = fit(x, k, fk, opts);

%% Visualize data

% Diffusion coefficients
figure;
subplot(3,1,1);
set(gcf,'color','w');
set(gca, 'YScale', 'log')
hold on;
plot(x(~exD), D(~exD), 'ko'); % Base points
plot(x(exD), D(exD), 'k+'); % Excluded points
plot(anode_fitD); % Fitting function
grid on;
legend({'Base points','Excluded points', 'Fit'},'Location','southwest');
xlabel('SOC') 
ylabel('Diffusion coefficient in m^2/s') 

% Anode potential
subplot(3,1,2);
plot(GITT.Results_Anode_20.OCV_q_rel_ch, GITT.Results_Anode_20.OCV_U_ch);
grid on;
xlabel('SOC') 
ylabel('Anode potential in V') 

% Reaction rates
subplot(3,1,3);
hold on;
plot(x(~exk), k(~exk), 'ko'); % Base points
plot(x(exk), k(exk), 'k+'); % Excluded points
plot(anode_fitk); % Fitting function
grid on;
legend({'Base points','Excluded points', 'Fit'},'Location','northwest');
xlabel('SOC') 
ylabel('Reaction rate in A m^{2.5}/mol^{1.5}')

%% Save data

savefig('final_results_kinetics');
save kinetics_anode anode_fitD anode_fitk