%% Curve fitting function of cathode OCV
%   Extended by Nikolaos Wassiliadis and Matthias Wanzel in December, 2019
%   Cathode fitting equation (polynomial function of 14th degree) to improve fitting quality

addpath(genpath(pwd));
clear;
close all;
clc;

%% Initialization
GITT = load('Halfcell_GITT'); % Load measurement
x = GITT.Results_Cathode_20.OCV_q_rel_ch; % Assign x-values
y = GITT.Results_Cathode_20.OCV_U_ch; % Assign y-values

%% Fitting function with 14th degree polynomial function
p = polyfit(x, y, 14); 

%% Calculate spline interpolation and regression for fit evaluation
x_ref = 0:1/10000:1; % Increase resolution for x
y_ref = interp1(x,y,x_ref,'spline'); % Spline interpolation for y

x_fit = x_ref; % Assign x
y_fit = polyval(p, x_fit); % Compute y regression fit

%% Calculate error of fit
r2 = 1 - sum((y_ref - y_fit).^2)/sum((y_ref - mean(y_ref)).^2); % R squared
NRMSE = goodnessOfFit(y_fit', y_ref','NRMSE'); % NRMSE

%%  Plot result
figure;
hold on;
plot(x,y,'b+'); % Measurement
plot(x_ref, y_ref, 'b'); % Spline interpolation of measurements
plot(x_fit, y_fit,'r'); % Plot regression fit
legend({'Reference','Reference interpolation', 'Regression'},'Location','best');

%% Save data
disp('Fitted function coefficients:')
disp(['p1 = ' num2str(p(1))]);
disp(['p2 = ' num2str(p(2))]);
disp(['p3 = ' num2str(p(3))]);
disp(['p4 = ' num2str(p(4))]);
disp(['p5 = ' num2str(p(5))]);
disp(['p6 = ' num2str(p(6))]);
disp(['p7 = ' num2str(p(7))]);
disp(['p8 = ' num2str(p(8))]);
disp(['p9 = ' num2str(p(9))]);
disp(['p10 = ' num2str(p(10))]);
disp(['p11 = ' num2str(p(11))]);
disp(['p12 = ' num2str(p(12))]);
disp(['p13 = ' num2str(p(13))]);
disp(['p14 = ' num2str(p(14))]);
disp(['p15 = ' num2str(p(15))]);
disp(' ');
disp(['R-squared yields: ' num2str(r2)])
disp(' ');

save OCV_cathode p