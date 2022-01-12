%% Create regression of half cell current step response fitting
%  Developed by Nikolaos Wassiliadis, Andreas Bach and Matthias Wanzel in November, 2020

addpath('./measurement');
addpath('./parameter');
addpath('./auxiliaries');
clear;
clc;
close all;

%% Load fitting (reaction rate and diffusion coefficient)
load('intermediate_results_coefficients');
GITT = load('Halfcell_GITT');

%% Fit regression

% Prepare curve fitting
[~,D] = prepareCurveData(coeff(:,1),coeff(:,2)); % Diffusion coefficients
[x,k] = prepareCurveData(coeff(:,1),coeff(:,3)); % Reaction rate

% Exclude errorenous data points, if necessary
exD = excludedata(x, D, 'Indices', [33 35:38]);
exk = excludedata(x, k, 'Indices', [33 35:38]);

% Define fitting functions
fD = fittype('poly3'); % Diffusion coefficient
fk = fittype('poly2'); % Reaction rate

% Fit
opts = fitoptions('Method', 'NonlinearLeastSquares', 'Algorithm', 'levenberg-marquardt');
opts.Exclude = exD; % Excluded points if necessary
[cathode_fitD, ~] = fit(x, D, fD, opts);
opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Algorithm', 'levenberg-marquardt');
opts.Exclude = exk; % Excluded points if necessary
[cathode_fitk, ~] = fit(x, k, fk, opts);

%% Visualize data

% Diffusion coefficients
figure;
subplot(2,1,1);
set(gcf,'color','w');
set(gca, 'YScale', 'log')
hold on;
plot(x(~exD), D(~exD), 'ko'); % Base points
plot(x(exD), D(exD), 'k+'); % Excluded points
plot(cathode_fitD); % Fitting function
grid on;
legend({'Base points','Excluded points', 'Fit'},'Location','southwest');
xlabel('SOC') 
ylabel('Diffusion coefficient in m^2/s')

% Reaction rates
subplot(2,1,2);
hold on; 
plot(x(~exk), k(~exk), 'ko'); % Base points
plot(x(exk), k(exk), 'k+'); % Excluded points
plot(cathode_fitk); % Fitting function
grid on;
legend({'Base points','Excluded points', 'Fit'},'Location','northwest');
xlabel('SOC') 
ylabel('Reaction rate in A m^{2.5}/mol^{1.5}') 

%% Save data
savefig('final_results_kinetics')
save kinetics_cathode cathode_fitD cathode_fitk