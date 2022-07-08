%% Sensitivity analysis of parameter variations to anode potential reserve control
%  Published in November, 2021 by Nikolaos Wassiliadis

addpath('./measurement');
addpath('./parameter');
addpath('./auxiliaries');
clear;
close all;
clc;
     
%% Load Datasets
SOH = 0.9; % Add SOH variation
SOC_start = 0.1;
SOC_end = 0.8;

% Get actual cell capacity
data.Q    = 2.5*SOH;

% Ensure sample frequency of 1Hz
t   = 0:0.1:(data.Q*(SOC_end-SOC_start))/15*3600; % 15A
I   = ones(1, length(t))*15;
SOC = (I.*t/3600/data.Q+SOC_start)*100;
T2  = ones(1, length(t))*20;

%% Configuration for sensitivity analysis

pscale = ones(16,1); % 
% pscale are the scaling factors used for parameter variation in this analysis
% All investigated parameter are grouped according to the occurence during parameter determination
% Group 1: Material and cell teardown
% pscale(1) - Active area
% pscale(2) - Particle radius of the anode
% pscale(3) - Particle radius of the cathode
% pscale(4) - Active material capacity of the anode
% pscale(5) - Active material capacity of the cathode
% Group 2: Open-circuit voltage determination
% pscale(6) - Shift (sigma) of anode open-circuit potential
% pscale(7) - Shift (sigma) of cathode open-circuit potential
% pscale(8) - Scaling (S) of anode open-circuit potential
% pscale(9) - Scaling (S) of cathode open-circuit potential
% Group 3: Dynamic parameter determination
% pscale(5) - Activation energy of the anode reaction rate
% pscale(6) - Activation energy of the cathode reaction rate
% pscale(7) - Activation energy of the anode diffusion coefficient
% pscale(8) - Activation energy of the cathode diffusion coefficient
% Group 4: Thermal parameter
% pscale(14) - Cell mass
% pscale(15) - Thermal capacity
% pscale(16) - Thermal resistance


%% Execute model

% Simulate new parameterization
[t_sim, U_sim_a, OCV_a] = spmet_fnc_halfcells(t, I, T2, data.Q, 20, 10, 5, 10, pscale);


%% Visualize: Anode potential over SOC
figure;
set(gcf,'color','w');
hold on;
grid on;
plot(SOC', U_sim_a);
xlabel('SOC in %')
ylabel('Anode potential in V')
xlim([10 80]);
ylim([0 0.1]);

% Downsampled for article publication
close all;
figure;
set(gcf,'color','w');
hold on;
grid on;
plot(SOC(1:10:end)', U_sim_a(1:10:end));
xlabel('SOC in %')
ylabel('Anode potential in V')
xlim([10 80]);
ylim([0 0.1]);
matlab2tikz('potential.tex')

% Downsampled for article publication
close all;
figure;
set(gcf,'color','w');
hold on;
grid on;
plot(SOC(1:10:end)', OCV_a(1:10:end));
xlabel('SOC in %')
ylabel('Anode potential in V')
xlim([10 80]);
ylim([0 0.1]);
matlab2tikz('potential_OCV.tex')


%% Execute model

% Parameter variation reach and step size
step = .01;
reach = .1;
window = (1-reach):step:(1+reach);
thres = 0.05; % 50mV

% Loop for data generation
for i=1:length(pscale)
    disp(['Running ' num2str(i) 'th line of parameter vector'])
    pscale = ones(16,1);
    for j=1:length(window)
        disp(['Sample ' num2str(j) ' out of ' num2str(length(window))])
        % Change pscale
        pscale(i)=window(j);
        % Simulate new parameterization
        [t_sim, U_sim_a] = spmet_fnc_halfcells(t, I, T2, data.Q, 20, 10, 5, 10, pscale);
        % Find undershoot
        try
        ix = find(U_sim_a<thres, 1, 'first');
        x = [U_sim_a(ix-1); U_sim_a(ix)];
        y = [SOC(ix-1); SOC(ix)];
        xq = thres;
        SOC_ix = interp1(x,y,xq);
        catch
        xq = NaN;
        SOC_ix = NaN;
        end
        % Transfer results
        results_simulation{i,j} = [SOC', U_sim_a];
        results_threshold{i,j} = [window(j), SOC_ix];
    end
end


%% Visualize: Sensitivity analysis

% Anode potential over SOC during variation
% figure;
% set(gcf,'color','w');
% hold on;
% for i=1:size(results_threshold,1)
%     for j=1:size(results_threshold,2)
%         subplot(size(results_simulation,1),1,i);
%         temp = cell2mat(results_simulation(i,j)');
%         plot(temp(:,1), temp(:,2),'+');
%         grid on;
%         xlabel('SOC in %')
%         ylabel('Potential in V')
%     end
% end

% Sensitivity of anode potential control point
figure;
set(gcf,'color','w');
hold on;
for i=1:size(results_threshold,1)
    % Realign results for visualization
    temp = cell2mat(results_threshold(i,:)');
    ix = find(temp(:,1)==1, 1, 'first');
    temp(:,1) = (temp(:,1)./temp(ix,1)-1)*100;
    temp(:,2) = (temp(:,2)./temp(ix,2)-1)*100;
    plot(temp(:,1),temp(:,2));
end
grid on;
xlabel('Parameter change in %')
ylabel('Charging limit position change in %')
legend('Area',...
       'Particle radius (anode)',...
       'Particle radius (cathode)',...       
       'Active material capacity (anode)',...
       'Active material capacity (cathode)',...        
       '\sigma (anode)',...
       '\sigma (cathode)',...
       'S (anode)',...
       'S (cathode)',...
       'E_k (anode)',...
       'E_k (cathode)',...
       'E_D (anode)',...
       'E_D (cathode)',...
       'Cell mass',...
       'Thermal capacity',...
       'Thermal resistance');

save('./results/analysis.mat');
savefig('./results/sensitivity.fig');


% Downsampled for article publication
close all;
figure;
set(gcf,'color','w');
hold on;
for i=1:size(results_threshold,1)
    % Realign results for visualization
    temp = cell2mat(results_threshold(i,:)');
    ix = find(temp(:,1)==1, 1, 'first');
    temp(:,1) = (temp(:,1)./temp(ix,1)-1)*100;
    temp(:,2) = (temp(:,2)./temp(ix,2)-1)*100;
    plot(temp(:,1),temp(:,2));
end
grid on;
xlabel('Parameter change in %')
ylabel('Charging limit position change in %')
legend('Area',...
       'Particle radius (anode)',...
       'Particle radius (cathode)',...       
       'Active material capacity (anode)',...
       'Active material capacity (cathode)',...        
       '\sigma (anode)',...
       '\sigma (cathode)',...
       'S (anode)',...
       'S (cathode)',...
       'E_k (anode)',...
       'E_k (cathode)',...
       'E_D (anode)',...
       'E_D (cathode)',...
       'Cell mass',...
       'Thermal capacity',...
       'Thermal resistance');
matlab2tikz('sensitivity.tex')