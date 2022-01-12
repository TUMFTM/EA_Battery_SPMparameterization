%% Objective function
%  See spmet.m for credentials

function [U] = spmet_objective(t, I, p0, SOC0, Nr, Nxn, Nxs, Nxp)

% Simulation of SPMeT
[~, ~, ~, ...
     ~, ~, ...
     ~, ~, ~, ...
     eta_n, ~, Unref, ~, ...
     eta_sei_n, ~, ...
     ~, ~] = spmet_fnc(t, I, p0, SOC0, Nr, Nxn, Nxs, Nxp);

% Return anode potential
U = real(Unref+eta_n+eta_sei_n);
end