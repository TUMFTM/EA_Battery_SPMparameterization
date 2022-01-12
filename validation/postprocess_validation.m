%% Evaluate model quality at different temperatures and C-Rates
%   Published in January, 2021 by Nikolaos Wassiliadis and Andreas Bach

clear;
close all;
clc;

%% Color definition - TUM colors

n = 50;                %// number of colors
TUM_map(1,:) = [255 255 255]./255;   %// color first row - TUM white
TUM_map(2,:) = [0 101 189]./255;   %// color 25th row - TUM blue
TUM_map(3,:) = [153 153 153]./255;   %// color 50th row - TUM gray

[X,Y] = meshgrid(1:3, 1:n);  %// mesh of indices

TUM_map = interp2(X([1,n/2,n],:),Y([1,n/2,n],:),TUM_map,X,Y); %// interpolate colormap


%% User input
% RMSE minimum of 10mV, RMSE maximum of 60mV
levels = 10:2:60;

%% Initialization

% Load data
temp     = load('results_RMSE_evaluation');
RMSE_new = temp.RMSE_eval;

temperatures = [0, 5, 10, 15, 20, 35, 50];
crates       = [1, 2, 3, 4, 5, 6];

%% Visualization

% New parameterization
figure('units','normalized','outerposition',[0 0 1 1]); % Full-screen figure with fitting results
set(gcf,'color','w');

contourf(temperatures, crates, RMSE_new, levels, 'ShowText', 'on');
h = colorbar;
ylabel(h,'RMSE in mV');
colormap(TUM_map)

yticks([1 2 3 4 5 6])
xlabel('Temperature in °C')
ylabel('C-Rate in h^{-1}')
title('Validation result - model quality')
savefig('results_model_quality')
