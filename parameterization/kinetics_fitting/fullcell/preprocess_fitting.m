%% Prepare dataset for pulse sequence fitting of full cell current step response
%  Developed by Nikolaos Wassiliadis and Andreas Bach in November, 2020

addpath('./auxiliaries');
addpath('./measurement');
addpath('./parameter');
clear;
close all;
clc;

%% Choose corresponding measurement data:

temperatures = [-10, 5, 20, 35, 50] + 273.15;

for idx_temp = 1:length(temperatures)

    % Define measurement temperature in °C and convert to K
    temp_measure  = temperatures(idx_temp);
    t_step        = 1;
    max_step      = 10;

    % Load GITT measurement
    fullcell = load(['measurement/fullcell_CP_' num2str(round(temp_measure,0)) '_Kelvin']);
    data     = fullcell.results_fullcell;

    % Load balancing and alignment
    B_A           = load('parameter/B_A.mat', 'Results');
    B_A           = B_A.Results;

    % Preallocate variable
    P = struct('t', cell(1, length(data.OCV_U_ch)-1), ...
               'U', cell(1, length(data.OCV_U_ch)-1), ...
               'I', cell(1, length(data.OCV_U_ch)-1));

    % Iterate through pulses and isolate to separate pulse files
    for i=1:length(data.OCV_U_ch)-1

        % Get time frame
        t_s = data.OCV_t_ch(i);
        t_e = data.OCV_t_ch(i+1);

        % Get indices of time frame
        ix_s = find(data.t_ch==t_s, 1, 'first');
        ix_e = find(data.t_ch<t_e, 1, 'last');

        % Get data
        t = data.t_ch(ix_s:ix_e-1);
        U = data.U_ch(ix_s:ix_e-1);
        I = data.I_ch(ix_s:ix_e-1);

        % Realign and resample to seconds
        [t, idx_t] = unique(t);
        t     = t-t(1);
        U     = U(idx_t);
        I     = I(idx_t);
        t_    = 0:1:t(end);
        U_uni = interp1(t, U, t_);
        I_uni = interp1(t, I, t_);

        % Assign to data
        P(i).t = t_;
        P(i).U = U_uni;
        P(i).I = I_uni;    
    end

    % Save results
    save(['intermediate_results_pulselist_' num2str(round(temp_measure,0)) '_Kelvin'], 'P', 'data')
end