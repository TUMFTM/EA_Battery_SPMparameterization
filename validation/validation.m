%% Validation of the parameterization technique
%   Published in January, 2021 by Nikolaos Wassiliadis and Andreas Bach

addpath('./measurement');
addpath('./parameter');
addpath('./auxiliaries');
clear;
close all;
clc;
     
%% Load Datasets

% Specify validation temperatures
temperatures = [-10, 0, 5, 10, 15, 20, 35, 50];

RMSE_eval     = zeros(6,8);

for idx_temp = 1:length(temperatures)

    temp_measure = temperatures(idx_temp);
    
    disp(['Constant Current validation at ' num2str(temp_measure) '°C'])

    % Specify number of measurements used for validation
    num_measure  = 6;

    % Preallocate cell
    measure = cell(num_measure, 1);

    % Get validation data
    if idx_temp == 1
        CC_1C        = load(['measurement/' num2str(temp_measure) '°C/US18650VTC5A_1C_m' num2str(abs(temp_measure)) 'deg']);
        measure{1,1} = CC_1C.Dataset;
        CC_2C        = load(['measurement/' num2str(temp_measure) '°C/US18650VTC5A_2C_m' num2str(abs(temp_measure)) 'deg']);
        measure{2,1} = CC_2C.Dataset;
        CC_3C        = load(['measurement/' num2str(temp_measure) '°C/US18650VTC5A_3C_m' num2str(abs(temp_measure)) 'deg']);
        measure{3,1} = CC_3C.Dataset;
        CC_4C        = load(['measurement/' num2str(temp_measure) '°C/US18650VTC5A_4C_m' num2str(abs(temp_measure)) 'deg']);
        measure{4,1} = CC_4C.Dataset;
        CC_5C        = load(['measurement/' num2str(temp_measure) '°C/US18650VTC5A_5C_m' num2str(abs(temp_measure)) 'deg']);
        measure{5,1} = CC_5C.Dataset;
        CC_6C        = load(['measurement/' num2str(temp_measure) '°C/US18650VTC5A_6C_m' num2str(abs(temp_measure)) 'deg']);
        measure{6,1} = CC_6C.Dataset;
    else
        CC_1C        = load(['measurement/' num2str(temp_measure) '°C/US18650VTC5A_1C_' num2str(temp_measure) 'deg']);
        measure{1,1} = CC_1C.Dataset;
        CC_2C        = load(['measurement/' num2str(temp_measure) '°C/US18650VTC5A_2C_' num2str(temp_measure) 'deg']);
        measure{2,1} = CC_2C.Dataset;
        CC_3C        = load(['measurement/' num2str(temp_measure) '°C/US18650VTC5A_3C_' num2str(temp_measure) 'deg']);
        measure{3,1} = CC_3C.Dataset;
        CC_4C        = load(['measurement/' num2str(temp_measure) '°C/US18650VTC5A_4C_' num2str(temp_measure) 'deg']);
        measure{4,1} = CC_4C.Dataset;
        CC_5C        = load(['measurement/' num2str(temp_measure) '°C/US18650VTC5A_5C_' num2str(temp_measure) 'deg']);
        measure{5,1} = CC_5C.Dataset;
        CC_6C        = load(['measurement/' num2str(temp_measure) '°C/US18650VTC5A_6C_' num2str(temp_measure) 'deg']);
        measure{6,1} = CC_6C.Dataset;
    end

    %% Merge datasets

    % Preallocate struct
    data = struct('t', cell(1, num_measure), ...
                  'U', cell(1, num_measure), ...
                  'I', cell(1, num_measure), ...
                  'Temp', cell(1, num_measure), ...
                  'Q', cell(1, num_measure));

    % Structure data in the right format
    for idx_data = 1:num_measure

        % Get actual cell capacity
        data(idx_data).Q    = abs(min(measure{idx_data,1}.AhStep));

        % Get relevant measurement range
        tmp = find(strcmp(measure{idx_data,1}.Command, 'CalcOnce')==1, 1, 'first');
        t   = measure{idx_data,1}.Time(tmp:end) - measure{idx_data,1}.Time(tmp);
        U   = measure{idx_data,1}.U(tmp:end);
        I   = measure{idx_data,1}.I(tmp:end);
        T1  = measure{idx_data,1}.T1(tmp:end);

        tmp = find(I~=0, 1, 'first')-1;
        t   = t(tmp:end) - t(tmp);
        U   = U(tmp:end);
        I   = I(tmp:end);
        T1  = T1(tmp:end);

        tmp = find(I~=0, 1, 'last');
        t   = t(1:tmp);
        U   = U(1:tmp);
        I   = I(1:tmp);
        T1  = T1(1:tmp);

        % Ensure sample frequency of 1Hz
        [t_uni, idx_uni] = unique(t*3600);
        tmp = 0:1:max(t_uni);
        U   = interp1(t_uni, U(idx_uni), tmp);
        I   = interp1(t_uni, I(idx_uni), tmp);
        T2  = interp1(t_uni, T1(idx_uni), tmp);

        data(idx_data).t    = tmp;
        data(idx_data).U    = U;
        data(idx_data).I    = I;
        data(idx_data).Temp = T2;
    end

    figure('units','normalized','outerposition',[0 0 1 1]); % Full-screen figure with fitting results
    set(gcf,'color','w');
    sgtitle(['Validation at ' num2str(temp_measure) '°C'])

    %% Simulate model for both parameterizations and every measurement

    rmse     = zeros(length(data),1);

    % Execute for every measurement loaded
    for idx_sim = 1:length(data)

        % Simulate new parameterization
        [t_sim, U_sim] = spmet_fnc(data(idx_sim).t, data(idx_sim).I, data(idx_sim).Temp, data(idx_sim).Q, ...
                                   20, 10, 5, 10);

        % Compute RMSE values in mV
        rmse(idx_sim,1)     = sqrt(mean((data(idx_sim).U'-U_sim).^2))*1000;
        
        RMSE_eval(idx_sim, idx_temp)     = sqrt(mean((data(idx_sim).U'-U_sim).^2))*1000;

        subplot(2,3,idx_sim)
        hold on
        grid on
        plot(t_sim, data(idx_sim).U, 'LineWidth', 1, 'Color', [0 0 0])
        plot(t_sim, U_sim, 'LineWidth', 1, 'Color', [162/255 173/255 0])
        xlabel('Time in s')
        ylabel('Voltage in V')
        legend('measurement', ['new fit: RMSE = ' num2str(round(rmse(idx_sim,1),2)) 'mV'], ...
               'Location', 'Southeast')

        switch num2str(idx_sim)
            case '1'
                title('Constant Current 1C')
            case '2'
                title('Constant Current 2C')
            case '3'
                title('Constant Current 3C')
            case '4'
                title('Constant Current 4C')
            case '5'
                title('Constant Current 5C')
            case '6'
                title('Constant Current 6C')
        end
    end
    
    if idx_temp == 1
        savefig(['results_validation_m' num2str(abs(temp_measure)) '°C'])
    else
        savefig(['results_validation_' num2str(temp_measure) '°C'])
    end
    
end

save('results_RMSE_evaluation', 'RMSE_eval')