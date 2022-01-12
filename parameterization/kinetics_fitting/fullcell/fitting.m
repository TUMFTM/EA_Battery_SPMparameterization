%% Pulse sequence fitting routine of full cell current step voltage response
%  Developed by Nikolaos Wassiliadis and Andreas Bach in November, 2020

addpath('./measurement');
addpath('./parameter');
addpath('./auxiliaries');
clear;
close all;
clc;
t1 = clock;

%% Iterate through all temperatures

temperatures = [-10, 5, 20, 35, 50] + 273.15; % define measurement temperatures in °C and convert to K

for idx_temp = 1:length(temperatures)
       
    % Initialization:
    % Adjust and assign input data
    temp_measure = temperatures(idx_temp); % Define measurement temperature
    
    disp(['Fitting at ' num2str(temp_measure-273.15) '°C'])

    % Define pulses to be excluded:
    switch num2str(round(temp_measure,0))
        case '263'
            idx_excluded = [4, 8, 12, 16, 20, 24, 28];
        case '278'
            idx_excluded = [4, 8, 12, 16, 20, 24, 28];
        case '293'
            idx_excluded = [4, 8, 12, 16, 20, 24, 28];
        case '308'
            idx_excluded = [4, 8, 12, 16, 20, 24, 28];
        case '323'
            idx_excluded = [4, 8, 12, 16, 20, 24, 28];
    end

    t_step       = 1; % Step size for SPM solving
    max_step     = 10; % 'MaxStep' option for ode23s solver

    % Load data 
    load(['intermediate_results_pulselist_' num2str(round(temp_measure,0)) '_Kelvin.mat']); % Load pulse list
    B_A    = load('parameter/B_A.mat', 'Results');
    B_A    = B_A.Results;

    %% Set solver specifications:

    % initial values used in the optimization algorithm lsqnonlin
    x0      = [1, 1, ...% [anodeK, cathodeK ...
               1, 1];  % [anodeD, cathodeD];

    % Print starting message  
    disp('Fitting started')

    % Prepare optimization
    options = optimoptions(@lsqnonlin,...
                           'Algorithm','levenberg-marquardt',...
                           'FiniteDifferenceType','forward',...
                           'FiniteDifferenceStepSize',1e-3,...
                           'TolX',  1e-25, ...
                           'Display', 'iter', ...
                           'TolFun', 1e-6, ...
                           'MaxIterations',1000, ...
                           'MaxFunEvals',30000,...
                           'UseParallel',true);

    %% Compute fitting

    % Define objective function:
    fun = @(x) spmet_objective(P, B_A, x, t_step, temp_measure, max_step, idx_excluded);

    % Start optimization:
    [x, resnorm, residual, exitflag, output] = lsqnonlin(fun, x0, [], [], options);

    % Save results
    x_opt = x;

    save(['intermediate_results_factors_' num2str(round(temp_measure,0)) '_Kelvin'], 'x_opt')

    %% Postprocessing

    % Show computation time
    t2 = clock;
    sim_time = etime(t2,t1);
    fprintf(1,'\n Parameter Fitting time: %4.1f sec or %2.2f min \n', sim_time, sim_time/60);

    % Visualize result
    figure('units','normalized','outerposition',[0 0 1 1]); % Full-screen figure with fitting results
    set(gcf,'color','w');

    pulses_RMSE_pre = zeros(length(P),1);
    pulses_RMSE_post = zeros(length(P),1); 

    for idx_pulse = 1:length(P)

        % Exclude long pulses used to change SOC levels
        if max(idx_pulse == idx_excluded)

            pulses_RMSE_pre(idx_pulse,1) = 0;
            pulses_RMSE_post(idx_pulse,1) = 0;

        else
            % Crawl data:
            t    = P(idx_pulse).t(1:t_step:end); % Downscale data to step size
            U    = P(idx_pulse).U(1:t_step:end); % Downscale data to step size
            I    = P(idx_pulse).I(1:t_step:end); % Downscale data to step size
            I(1) = eps; % Reset initial current to zero for ODE-solver compability

            % Limit data to samples where I~=0 plus same length of relaxation
            temp = find(I>0, 1, 'last');
            temp = 2*temp;
            t    = t(1:temp);
            U    = U(1:temp);
            I    = I(1:temp);

            % Find initial SOC
            U0 = U(1); % Hand-over initial OCV
            [OCV_uni, ix] = unique(B_A.OCV_fullcell_fit); % Prep data for inverse SOC-OCV look-up
            %
            if length(ix) > length(B_A.q_rel_dot_fullcell)
                ix      = ix(1:length(B_A.q_rel_dot_fullcell));
                OCV_uni = OCV_uni(1:length(B_A.q_rel_dot_fullcell));
            end
            %
            SOC0 = interp1(OCV_uni, B_A.q_rel_dot_fullcell(ix), U0, 'spline'); % Find SOC with inverse SOC-OCV look-up

            % Simulate pulse
            [U_sim_pre] = spmet_fnc(x0, t, I, temp_measure, SOC0, 20, 10, 5, 10, max_step);
            [U_sim_post] = spmet_fnc(x, t, I, temp_measure, SOC0, 20, 10, 5, 10, max_step);

            % Compute root mean square error in mV
            pulses_RMSE_pre(idx_pulse,1) = sqrt(mean((U'-U_sim_pre).^2));
            pulses_RMSE_post(idx_pulse,1)= sqrt(mean((U'-U_sim_post).^2));

            subplot(4,7,idx_pulse)
            hold on
            grid on
            title(['Pulse ' num2str(idx_pulse)])
            plot(t, U, 'LineWidth', 1, 'Color', [0 0 0])
            plot(t, U_sim_pre, 'LineWidth', 1, 'Color', [0/255 101/255 189/255])
            plot(t, U_sim_post,'LineWidth', 1, 'Color', [162/255 173/255 0/255])

        end

    end

    savefig(['intermediate_results_pulses_' num2str(round(temp_measure,0)) '_Kelvin'])
end