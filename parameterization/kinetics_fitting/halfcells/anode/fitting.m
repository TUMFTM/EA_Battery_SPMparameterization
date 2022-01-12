%% Pulse fitting routine of half cell current step voltage responses
%  Developed by Nikolaos Wassiliadis, Andreas Bach and Matthias Wanzel in November, 2020

addpath('./measurement');
addpath('./parameter');
addpath('./auxiliaries');
clear;
clc;
close all;

%% Initialize

% Load data
load('intermediate_results_pulselist.mat'); % Load pulse list
B_A = load('B_A.mat'); % Load balancing and alignment parameter for half cell scaling
B_A = B_A.Results;
t_step = 10; % Step size for SPM solving

% Visualize results
figure('units','normalized','outerposition',[0 0 1 1]); % Full-screen figure with fitting results
set(gcf,'color','w');

%% Start fitting

coeff = zeros(size(P,2)-1,3); % Preallocate results

for i=1:size(P,2)-1
    disp(['Pulse ' num2str(i) ':']);
    tic;

    % Crawl data
    I = P(i).I(1:t_step:end)*-1*B_A.s_neg; % Inverse and scale for anode in full cell simulation setup    
    t = P(i).t(1:t_step:end); % Downscale data to step size
    U = P(i).U(1:t_step:end); % Downsacle data to step size
    temp = find(I~=0, 1, 'last') * 2; % Limit data to samples where I~=0 plus same time of relaxation
    I = I(1:temp); % Downscale data to step size
    t = t(1:temp); % Downscale data to step size
    U = U(1:temp); % Downsacle data to step size
    I(1) = eps; % Reset initial current to zero for ODE-solver compability          
    
    % Find initial SOC for precise determination of OCV
    SOC0 = P(i).SOC;
        
    % Prep fitting
    p0 = [1, ... % factor to manipulate diffusion coeff for solid in neg. electrode
          1];    % factor to manipulate reaction rate in neg. electrode
    
    options = optimoptions(@lsqnonlin,...
                       'Algorithm','levenberg-marquardt',... % Bound-less fast solver
                       'FiniteDifferenceStepSize', 1e-3,... % Increased finite step size to default
                       'TolX',  1e-25, ... % Minimum step tolerance to terminate fitting procedure
                       'Display', 'iter-detailed', ... % Detailed optimization protocol in command
                       'TolFun', 1e-6, ... % Minimum step delta objective function
                       'MaxIterations',30, ... % Limit iterations for better performance
                       'MaxFunEvals',500, ... % Limit function evaluations for better performance
                       'UseParallel',true); % Increase performance with parallel computing (optional)
    
    f = @(p) U' - spmet_objective(t', I', p, SOC0, 20, 10, 5, 10); % Set up objective function
    
    % Simulate pulse and visualize before fitting
    [t_sim, ~, ~, ...
     ~, ~, ...
     ~, ~, ~, ...
     eta_n, ~, Unref, ~, ...
     eta_sei_n, ~, ...
     ~, ~] = spmet_fnc(t, I, p0, SOC0, 20, 10, 5, 10);

    subplot(4,10,i);
    hold on;
    plot(t, U, 'k');
    plot(t, real(Unref+eta_n+eta_sei_n), 'Color', [152 198 234]./255);
    
    % Start fitting procedure
    [p, resnorm, residual, exitflag, output] = lsqnonlin(f, p0, [], [], options);
    coeff(i,:) = [SOC0 p];
    
    % Simulate pulse, calculate RMSE and visualize after fitting
    [~, ~, ~, ...
     SOC_n, SOC_p, ...
     ~, ~, ~, ...
     eta_n, ~, Unref, ~, ...
     eta_sei_n, ~, ...
     ~, ~] = spmet_fnc(t, I, p, SOC0, 20, 10, 5, 10);
 
    RMSE = sqrt(mean((U' - real(Unref+eta_n+eta_sei_n)).^2));  % Root Mean Squared Error
 
    subplot(4,10,i);
    hold on;
    plot(t, real(Unref+eta_n+eta_sei_n), 'Color', [227 114 34]./255);
    title(['RMSE = ' num2str(round(RMSE*1000,2)) 'mV'])
    toc;
end

%% Save results
savefig('intermediate_results_fitting');

coeff(:,2) = coeff(:,2) * 7e-15;
coeff(:,3) = coeff(:,3) * 3e-6;
save intermediate_results_coefficients coeff;