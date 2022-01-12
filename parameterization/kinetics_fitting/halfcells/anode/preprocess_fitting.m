%% Prepare dataset for pulse fitting of half cell current step responses
%  Developed by Nikolaos Wassiliadis, Andreas Bach and Matthias Wanzel in November, 2020

addpath('./measurement');
addpath('./parameter');
addpath('./auxiliaries');
clear;
clc;
close all;

%% Isolate pulses

% Get and assign measurement data
Halfcell_GITT = load('measurement/Halfcell_GITT'); % Load GITT measurement
data  = Halfcell_GITT.Results_Anode_20; % Use 20°C measurement

% Preallocate variable
P = struct('t', cell(1, length(data.OCV_U_ch)-1), ...
           'U', cell(1, length(data.OCV_U_ch)-1), ...
           'I', cell(1, length(data.OCV_U_ch)-1), ...
           'SOC', cell(1, length(data.OCV_U_ch)-1));

% Iterate through pulses and isolate to separate pulse files
for i=1:length(data.OCV_U_ch)-1
    
    % Get time frame
    t_s = data.OCV_t_ch(i);
    t_e = data.OCV_t_ch(i+1);
    
    % Get indices of time frame
    ix_s = find(data.t_ch==t_s,1, 'first');
    ix_e = find(data.t_ch<t_e,1,'last');
    
    % Get data
    t = data.t_ch(ix_s:ix_e-1)*3600; % Scale to seconds
    U = data.U_ch(ix_s:ix_e-1);
    I = data.I_ch(ix_s:ix_e-1);
        
    % Realign and resample to seconds
    t = t-t(1); 
    t_uni = 0:1:t(end);
    U_uni = interp1(t,U,t_uni);
    I_uni = interp1(t,I,t_uni);
    
    % Assign to data
    P(i).t = t_uni;
    P(i).U = U_uni;
    P(i).I = I_uni;    
    P(i).SOC = data.OCV_q_rel_ch(i);
end

% Save results
save intermediate_results_pulselist P data