%% Objective function of fullcell k & D fitting
%  Developed by Nikolaos Wassiliadis and Andreas Bach in January, 2021

function [sum_RMSE] = spmet_objective(P, B_A, x, t_step, temp_measure, max_step, idx_excluded)
% P = pulse list
% B_A = OCV regression functions for anode and cathode
% x0 = initial point used in lsqnonlin
% options = solver options for lsqnonlin
% t_step = time discretization
% temp_measure = temperature of measurement data

%% Preallocation and Iteration through all pulses
disp('New function evaluation')

pulses_RMSE   = zeros(length(P),1);

for idx_pulse = 1:length(P)
    
    if max(idx_pulse == idx_excluded)
   
        pulses_RMSE(idx_pulse,1)   = 0;
        
    else
        % Crawl data:
        t    = P(idx_pulse).t(1:t_step:end); % Downscale data to step size
        U    = P(idx_pulse).U(1:t_step:end); % Downscale data to step size
        I    = P(idx_pulse).I(1:t_step:end); % Downscale data to step size
        I(1) = 0; % Reset initial current to zero for ODE-solver compability

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
        [U_sim] = spmet_fnc(x, t, I, temp_measure, SOC0, 20, 10, 5, 10, max_step);
        
        pulses_RMSE(idx_pulse,1)   = sqrt(mean((U'-U_sim).^2));
        
    end
    
end

sum_RMSE = sum(pulses_RMSE);

end

