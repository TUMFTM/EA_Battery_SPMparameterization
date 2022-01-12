%% Sensitivity analysis of parameter variations to anode potential reserve control
%  Published in January, 2022 by Nikolaos Wassiliadis

addpath('./measurement');
addpath('./parameter');
addpath('./auxiliaries');
clear;
close all;
clc;
     
%% Load Datasets
load('measurement/US18650VTC5A_6C_0deg');

% Get actual cell capacity
data.Q    = abs(min(Dataset.AhStep));

% Get relevant measurement range
tmp = find(strcmp(Dataset.Command, 'CalcOnce')==1, 1, 'first');
t   = Dataset.Time(tmp:end) - Dataset.Time(tmp);
U   = Dataset.U(tmp:end);
I   = Dataset.I(tmp:end);
AhStep = Dataset.AhStep(tmp:end);
 
% Realign time
tmp = find(I~=0, 1, 'first')-1;
t   = t(tmp:end) - t(tmp);
U   = U(tmp:end);
I   = I(tmp:end);
AhStep = AhStep(tmp:end);
 
tmp = find(I~=0, 1, 'last');
t   = t(1:tmp);
U   = U(1:tmp);
I   = I(1:tmp);
AhStep = AhStep(1:tmp);
 
% Ensure sample frequency of 10Hz
[t_uni, idx_uni] = unique(t*3600);
tmp = 0:0.1:max(t_uni);
t   = tmp;
U   = interp1(t_uni, U(idx_uni), tmp);
I   = interp1(t_uni, I(idx_uni), tmp);
T2  = zeros(1,length(I));
SOC = interp1(t_uni, AhStep(idx_uni), tmp)*100./data.Q;

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

for idx_sim = 1:length(data)
% Simulate new parameterization
[t_sim, U_sim_a] = spmet_fnc_halfcells(t, I, T2, data.Q, 20, 10, 5, 10, pscale);
end


%% Visualize: Anode potential over SOC
figure;
set(gcf,'color','w');
hold on;
grid on;
plot(SOC', U_sim_a);
xlabel('SOC in %')
ylabel('Anode potential in V')
xlim([0 80]);
ylim([0 0.1]);

% Downsampled for article publication
%close all;
figure;
set(gcf,'color','w');
hold on;
grid on;
plot(SOC(1:10:end)', U_sim_a(1:10:end));
xlabel('SOC in %')
ylabel('Anode potential in V')
xlim([0 80]);
ylim([0 0.1]);
%matlab2tikz('potential.tex')

%% Execute model

% Parameter variation reach and step size
step = .1;
reach = .1;
window = (1-reach):step:(1+reach);
thres = 0.05; % 50mV

figure;
hold on;
% Loop for data generation
for i=1:1%length(pscale)
    pscale = ones(16,1);
    for j=1:length(window)
        % Change pscale
        pscale(i)=window(j);
        % Simulate new parameterization
        [t_sim, U_sim_a] = spmet_fnc_halfcells(t, I, T2, data.Q, 20, 10, 5, 10, pscale);
        % Find undershoot
        ix = find(U_sim_a<thres, 1, 'first');
        x = [U_sim_a(ix-1); U_sim_a(ix)];
        y = [SOC(ix-1); SOC(ix)];
        xq = thres;
        SOC_ix = interp1(x,y,xq);
        % Transfer results
        results_simulation{i,j} = [SOC', U_sim_a];
        results_threshold{i,j} = [window(j), SOC_ix];
        plot(SOC',U_sim_a);
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
       'Active material capacity (cathide)',...
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
%close all;
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
       'Active material capacity (cathide)',...
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