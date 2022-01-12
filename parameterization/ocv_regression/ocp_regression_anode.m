%% Curve fitting function of anode OCV
%   Extended by Nikolaos Wassiliadis and Matthias Wanzel in December, 2019
%   Anode fitting equation taken from Safari and Delacourt (2011) - DOI: 10.1149/1.3567007

addpath(genpath(pwd));
clear;
close all;
clc;

%% Initialization
GITT = load('Halfcell_GITT'); % Load measurement
x = GITT.Results_Anode_20.OCV_q_rel_ch; % Assign x-values
y = GITT.Results_Anode_20.OCV_U_ch; % Assign y-values

%% Fitting function
p0 = [ 0.6379,... % Starting points from Safari and Delacourt
       0.5416,...
       -305.5309,...
       0.044,...
       0.1958,...
       -0.1088,...
       -0.1978,...
       -1.0571,...
       0.0854,...
       -0.6875,...
       0.0117,...
       0.0529,...
       -0.0175,...
       -0.5692,...
       0.0875];

fun = @(p)(p(1)+p(2)*exp(p(3)*x)+p(4)*tanh((x+p(5))/p(6))+p(7)*tanh((x+p(8))/p(9))+p(10)*tanh((x+p(11))/p(12))+p(13)*tanh((x+p(14))/p(15))) - y; % Objective function
options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter','MaxFunctionEvaluations',25000,'MaxIterations',10000,'FunctionTolerance',1e-16,'StepTolerance', 1e-20); % Optimization parameter
[p,resnorm,residual,exitflag,output] = lsqnonlin(fun,p0,[],[],options); % Optimization

%% Calculate spline interpolation and regression for fit evaluation
x_ref = 0:1/10000:1; % Increase resolution for x
y_ref = interp1(x,y,x_ref,'spline'); % Spline interpolation for y

x_fit = x_ref; % Assign x
y_fit = p(1)+p(2)*exp(p(3)*x_fit)+p(4)*tanh((x_fit+p(5))/p(6))+p(7)*tanh((x_fit+p(8))/p(9))+p(10)*tanh((x_fit+p(11))/p(12))+p(13)*tanh((x_fit+p(14))/p(15)); % Compute y regression fit

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

save OCV_anode p


