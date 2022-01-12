%% Arrhenius fitting
%   Developed by Nikolaos Wassiliadis and Andreas Bach in November, 2020.

function [coeff] = arrh_createFit(x, y, parameter)

%CREATEFIT(X_K_N,Y_K_N)
%  Create a fit.
%
%  Data for fit:
%      x        : input xData
%      y        : input yData
%      parameter: parameter specification
%      
%  Output:
%      coeff    : coefficients at reference temperature 293.15K 


%% Fitting:

% Set up fittype and options
fun     = @(coeff, x) coeff(1) + coeff(2)/8.314472 * (1/293.15-1./x);
x0      = [1 1];
options = optimoptions('lsqcurvefit', 'Algorithm', 'levenberg-marquardt');

% Fit model to data.
[coeff] = lsqcurvefit(fun, x0, x, log(y), [], [], options);

x_fit = 253.15:0.01:333.15;

fun_fit   = coeff(1) + coeff(2)/8.314472 * (1/293.15-1./x_fit);

coeff(1)  = exp(coeff(1));

%% Visualization

% Plot fit with data
switch parameter
    
    case 'k_n'
        figure('units','normalized','outerposition',[0 0 1 1]); % Full-screen figure with fitting results
        set(gcf,'color','w');
        subplot(2,2,1)
        
    case 'k_p'
        subplot(2,2,2)
        
    case 'D_n'
        subplot(2,2,3)
        
    case 'D_p'
        subplot(2,2,4)
        
end

plot(1./x, log(y), 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 6, 'Color', [0 0 0])
hold on
plot(1./x_fit, fun_fit, 'LineStyle', '-', 'LineWidth', 1.2, 'Color', [0 101/255 189/255])
grid on

% Label axes
xlabel('1/Temperature in 1/K');
ylabel('log(Scaling factor)');
title([parameter ' Arrhenius fitting'])
legend('data', 'fit', 'Location', 'northeast')

end
        