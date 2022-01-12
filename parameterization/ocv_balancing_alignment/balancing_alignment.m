%% Balancing and alignment of half-cells to a full-cell
%   Method is similar to the method published by Honkura and Horiba in https://doi.org/10.1016/j.jpowsour.2014.04.036
%   Published in November 2020 by Nikolaos Wassiliadis and Matthias Wanzel

addpath(genpath(pwd));

%% Selection of different measurement techniques
% Mode 1 - Halfcell GITT         + Fullcell GITT       with DVA objective function
% Mode 2 - Halfcell GITT fit     + Fullcell GITT       with DVA objective function
% Mode 3 - Halfcell GITT         + Fullcell GITT       with OCV objective function
% Mode 4 - Halfcell GITT fit     + Fullcell GITT       with OCV objective function < deployed
mode = 4;

%% Main routine with balancing and alignment of half cell potential curves to the full cell potential

% Load OCV regressions into workspace
transfer = load('OCV_anode'); % Anode GITT regression
p_a = transfer.p;
transfer = load('OCV_cathode'); % Cathode GITT regression
p_c = transfer.p;

% Load OCV measurements
Halfcell_GITT = load('Halfcell_GITT'); % Half cell GITT measurement
Fullcell_GITT = load('Fullcell_GITT'); % Full cell GITT measurement

% Configure balancing and alignment
reso = 10000; % Fitting resolution   
DVAstart = 0.05; % Considered DVA range starting point as SOC
DVAend = 0.80; % Considered DVA range ending point as SOC

% Mode 1 - Halfcell GITT -> fullcell GITT with DVA objective function
if mode == 1
    % Prepare full cell data to chosen resolution step size
    Q_fullcell = (0:max(Fullcell_GITT.Results_Fullcell_20.q_ch)/reso:max(Fullcell_GITT.Results_Fullcell_20.q_ch))*3600; % Scale full cell charge
    OCV_fullcell = interp1(Fullcell_GITT.Results_Fullcell_20.OCV_q_ch*3600, Fullcell_GITT.Results_Fullcell_20.OCV_U_ch, Q_fullcell,'spline'); % Scale full cell OCV(by spline interpolation)
    q_dot_fullcell = Q_fullcell(2:length(Q_fullcell)); % Calculate first derivative of full cell charge for DVA
    u_dot_fullcell = diff(OCV_fullcell)./diff(Q_fullcell); % Calculate first derivative of full cell OCV for DVA
    u_rel_dot_fullcell = diff(OCV_fullcell)./diff(Q_fullcell/max(Q_fullcell)); % Calculate normalized first derivative of full cell OCV for DVA
    q_rel_dot_fullcell = (Q_fullcell(2:length(Q_fullcell))/max(Q_fullcell)); % Calculate normalized first derivative of full cell charge for DVA

    % Prepare half cell data
    OCV_anode = Halfcell_GITT.Results_Anode_20.OCV_U_ch;
    OCV_cathode = Halfcell_GITT.Results_Cathode_20.OCV_U_ch;
    Q_anode = (0:max(abs(Halfcell_GITT.Results_Anode_20.OCV_q_ch))/reso:max(abs(Halfcell_GITT.Results_Anode_20.OCV_q_ch)))*3600; % Scale anode charge
    Q_cathode = (0:max(Halfcell_GITT.Results_Cathode_20.OCV_q_ch)/reso:max(Halfcell_GITT.Results_Cathode_20.OCV_q_ch))*3600; % Scale cathode charge
    [x, index] = unique(abs(Halfcell_GITT.Results_Anode_20.OCV_q_ch*3600));
    OCV_anode = interp1(x, OCV_anode(index), Q_anode,'spline');
    [x, index] = unique(Halfcell_GITT.Results_Cathode_20.OCV_q_ch*3600);
    OCV_cathode = interp1(x, OCV_cathode(index), Q_cathode,'spline');

    % Get relevant DVA range
    rg = find(DVAstart*max(Q_fullcell)<Q_fullcell & Q_fullcell<DVAend*max(Q_fullcell));
    
    % Balancing and alignment (objective function)
    % Optimization routine minimizes error between measured OCV and calculated OCV from both electrode potentials
    % x(1) and x(2) offset electrode capacities
    % x(3) and x(4) scale electrode capacities
    fun = @(x) ...
        (diff(OCV_fullcell(rg))./diff(Q_fullcell(rg))) ... % Full-cell DVA
        -(diff(interp1(Q_cathode*x(4)-x(2), OCV_cathode, Q_fullcell(rg),'linear') ... % Cathode DVA
        -interp1(Q_anode*x(3)-x(1), OCV_anode, Q_fullcell(rg),'linear'))./diff(Q_fullcell(rg))); % Anode DVA

% Mode 2 - Halfcell GITT fitted -> fullcell GITT with DVA objective function
elseif mode == 2 
    % Prepare full cell data to chosen resolution step size
    Q_fullcell = (0:max(Fullcell_GITT.Results_Fullcell_20.q_ch)/reso:max(Fullcell_GITT.Results_Fullcell_20.q_ch))*3600; % Scale full cell charge
    OCV_fullcell = interp1(Fullcell_GITT.Results_Fullcell_20.OCV_q_ch*3600, Fullcell_GITT.Results_Fullcell_20.OCV_U_ch, Q_fullcell,'spline'); % Scale full cell OCV(by spline interpolation)
    q_dot_fullcell = Q_fullcell(2:length(Q_fullcell)); % Calculate first derivative of full cell charge for DVA
    u_dot_fullcell = diff(OCV_fullcell)./diff(Q_fullcell); % Calculate first derivative of full cell OCV for DVA
    u_rel_dot_fullcell = diff(OCV_fullcell)./diff(Q_fullcell/max(Q_fullcell)); % Calculate normalized first derivative of full cell OCV for DVA
    q_rel_dot_fullcell = (Q_fullcell(2:length(Q_fullcell))/max(Q_fullcell)); % Calculate normalized first derivative of full cell charge for DVA
    
    % Prepare half cell data to chosen resolution step size
    theta = 0:1/reso:1; % Calculate normalized resolution
    OCV_anode = p_a(1)+p_a(2)*exp(p_a(3)*theta)+p_a(4)*tanh((theta+p_a(5))/p_a(6))+p_a(7)*tanh((theta+p_a(8))/p_a(9))+p_a(10)*tanh((theta+p_a(11))/p_a(12))+p_a(13)*tanh((theta+p_a(14))/p_a(15)); % Scale anode OCV (hyperbolic fit)
    OCV_cathode = polyval(p_c, theta); % Scale cathode OCV (poly-fit)
    Q_anode = (0:max(abs(Halfcell_GITT.Results_Anode_20.q_ch))/reso:max(abs(Halfcell_GITT.Results_Anode_20.q_ch)))*3600; % Scale anode charge
    Q_cathode = (0:max(Halfcell_GITT.Results_Cathode_20.q_ch)/reso:max(Halfcell_GITT.Results_Cathode_20.q_ch))*3600; % Scale cathode charge

    % Get relevant DVA range
    rg = find(DVAstart*max(Q_fullcell)<Q_fullcell & Q_fullcell<DVAend*max(Q_fullcell));
    
    % Balancing and alignment (objective function)
    % Optimization routine minimizes error between measured OCV and calculated OCV from both electrode potentials
    % x(1) and x(2) offset electrode capacities
    % x(3) and x(4) scale electrode capacities
    fun = @(x) ...
        (diff(OCV_fullcell(rg))./diff(Q_fullcell(rg))) ... % Full-cell DVA
        -(diff(interp1(Q_cathode*x(4)-x(2), OCV_cathode, Q_fullcell(rg),'linear') ... % Cathode DVA
        -interp1(Q_anode*x(3)-x(1), OCV_anode, Q_fullcell(rg),'linear'))./diff(Q_fullcell(rg))); % Anode DVA

% Mode 3 - Halfcell GITT -> fullcell GITT with DVA objective function
elseif mode == 3
    % Prepare full cell data to chosen resolution step size
    Q_fullcell = (0:max(Fullcell_GITT.Results_Fullcell_20.q_ch)/reso:max(Fullcell_GITT.Results_Fullcell_20.q_ch))*3600; % Scale full cell charge
    OCV_fullcell = interp1(Fullcell_GITT.Results_Fullcell_20.OCV_q_ch*3600, Fullcell_GITT.Results_Fullcell_20.OCV_U_ch, Q_fullcell,'spline'); % Scale full cell OCV(by spline interpolation)
    q_dot_fullcell = Q_fullcell(2:length(Q_fullcell)); % Calculate first derivative of full cell charge for DVA
    u_dot_fullcell = diff(OCV_fullcell)./diff(Q_fullcell); % Calculate first derivative of full cell OCV for DVA
    u_rel_dot_fullcell = diff(OCV_fullcell)./diff(Q_fullcell/max(Q_fullcell)); % Calculate normalized first derivative of full cell OCV for DVA
    q_rel_dot_fullcell = (Q_fullcell(2:length(Q_fullcell))/max(Q_fullcell)); % Calculate normalized first derivative of full cell charge for DVA
    
    % Prepare half cell data
    OCV_anode = Halfcell_GITT.Results_Anode_20.OCV_U_ch;
    OCV_cathode = Halfcell_GITT.Results_Cathode_20.OCV_U_ch;
    Q_anode = (0:max(abs(Halfcell_GITT.Results_Anode_20.OCV_q_ch))/reso:max(abs(Halfcell_GITT.Results_Anode_20.OCV_q_ch)))*3600; % Scale anode charge
    Q_cathode = (0:max(Halfcell_GITT.Results_Cathode_20.OCV_q_ch)/reso:max(Halfcell_GITT.Results_Cathode_20.OCV_q_ch))*3600; % Scale cathode charge
    [x, index] = unique(abs(Halfcell_GITT.Results_Anode_20.OCV_q_ch*3600));
    OCV_anode = interp1(x, OCV_anode(index), Q_anode,'spline');
    [x, index] = unique(Halfcell_GITT.Results_Cathode_20.OCV_q_ch*3600);
    OCV_cathode = interp1(x, OCV_cathode(index), Q_cathode,'spline');

    %Balancing and alignment (objective function)
    % Optimization routine minimizes error between measured OCV and calculated OCV from both electrode potentials
    % x(1) and x(2) offset electrode capacities
    % x(3) and x(4) scale electrode capacities
    fun = @(x) ...
        OCV_fullcell ... % Full-cell OCV
        -(interp1(Q_cathode*x(4)-x(2), OCV_cathode, Q_fullcell,'linear') ... % Cathode OCV
        -interp1(Q_anode*x(3)-x(1), OCV_anode, Q_fullcell,'linear')); % Anode OCV

% Mode 4 - Halfcell GITT fitted -> fullcell GITT with OCV objective function
elseif mode == 4
    % Prepare full cell data to chosen resolution step size
    Q_fullcell = (0:max(Fullcell_GITT.Results_Fullcell_20.q_ch)/reso:max(Fullcell_GITT.Results_Fullcell_20.q_ch))*3600; % Scale full cell charge
    OCV_fullcell = interp1(Fullcell_GITT.Results_Fullcell_20.OCV_q_ch*3600, Fullcell_GITT.Results_Fullcell_20.OCV_U_ch, Q_fullcell,'spline'); % Scale full cell OCV(by spline interpolation)
    q_dot_fullcell = Q_fullcell(2:length(Q_fullcell)); % Calculate first derivative of full cell charge for DVA
    u_dot_fullcell = diff(OCV_fullcell)./diff(Q_fullcell); % Calculate first derivative of full cell OCV for DVA
    u_rel_dot_fullcell = diff(OCV_fullcell)./diff(Q_fullcell/max(Q_fullcell)); % Calculate normalized first derivative of full cell OCV for DVA
    q_rel_dot_fullcell = (Q_fullcell(2:length(Q_fullcell))/max(Q_fullcell)); % Calculate normalized first derivative of full cell charge for DVA
    
    % Prepare half cell data to chosen resolution step size
    theta = 0:1/reso:1; % Calculate normalized resolution
    OCV_anode = p_a(1)+p_a(2)*exp(p_a(3)*theta)+p_a(4)*tanh((theta+p_a(5))/p_a(6))+p_a(7)*tanh((theta+p_a(8))/p_a(9))+p_a(10)*tanh((theta+p_a(11))/p_a(12))+p_a(13)*tanh((theta+p_a(14))/p_a(15)); % Scale anode OCV (hyperbolic fit)
    OCV_cathode = polyval(p_c, theta); % Scale cathode OCV (poly-fit)
    Q_anode = (0:max(abs(Halfcell_GITT.Results_Anode_20.q_ch))/reso:max(abs(Halfcell_GITT.Results_Anode_20.q_ch)))*3600; % Scale anode charge
    Q_cathode = (0:max(Halfcell_GITT.Results_Cathode_20.q_ch)/reso:max(Halfcell_GITT.Results_Cathode_20.q_ch))*3600; % Scale cathode charge

    % Balancing and alignment (objective function)
    % Optimization routine minimizes error between measured OCV and calculated OCV from both electrode potentials
    % x(1) and x(2) offset electrode capacities
    % x(3) and x(4) scale electrode capacities
    fun = @(x) ...
        OCV_fullcell ... % Full-cell OCV
        -(interp1(Q_cathode*x(4)-x(2), OCV_cathode, Q_fullcell,'linear') ... % Cathode OCV
        -interp1(Q_anode*x(3)-x(1), OCV_anode, Q_fullcell,'linear')); % Anode OCV
end

% Set starting parameter for optimization
sigma_neg = 2; % Negative shift
sigma_pos = 1250; % Positive shift
s_neg = 530; % Negative scale
s_pos = 505; % Positive scale
x0 = [sigma_neg, sigma_pos, s_neg, s_pos]; % Assign to optimization

% Optimization configuration
options = optimoptions(@lsqnonlin, ...
    'Algorithm','levenberg-marquardt', ...
    'Display','iter-detailed', ... 
    'FiniteDifferenceType','central', ...
    'MaxFunctionEvaluations',50000, ...
    'FunctionTolerance', 1e-18, ...
    'StepTolerance', 1e-24);

% Objective function hand-over
x = lsqnonlin(fun, x0,[],[] , options);

% Fitting parameters
Results.sigma_neg = x(1);
Results.sigma_pos = x(2);
Results.s_neg = x(3);
Results.s_pos = x(4);

% Capacity and OCV of fullcell 
Results.Q_fullcell = Q_fullcell;
Results.OCV_fullcell = OCV_fullcell;

% Capacity and OCV of fitted half cells
Results.Q_anode = Q_anode*x(3)-x(1);
Results.Q_cathode = Q_cathode*x(4)-x(2);
Results.OCV_anode = OCV_anode;
Results.OCV_cathode = OCV_cathode;

% Used half cell OCV range (q = Q_fullcell)
Results.OCV_anode_fit = interp1(Q_anode*x(3)-x(1), OCV_anode, Q_fullcell,'linear');
Results.OCV_cathode_fit = interp1(Q_cathode*x(4)-x(2), OCV_cathode, Q_fullcell,'linear');

% Superposition of the used half cell OCV range  (q = Q_fullcell)
Results.OCV_fullcell_fit = interp1(Q_cathode*x(4)-x(2), OCV_cathode, Q_fullcell,'linear')-interp1(Q_anode*x(3)-x(1), OCV_anode, Q_fullcell,'linear');

% Relative start- and endpoint of the used capacity ranges of the fitted half cell capacities
Results.Q_anode_SP = (min(Q_fullcell)-min(Q_anode*x(3)-x(1)))/(max(Q_anode*x(3)-x(1))-min(Q_anode*x(3)-x(1)));
Results.Q_anode_EP = (max(Q_fullcell)-min(Q_anode*x(3)-x(1)))/(max(Q_anode*x(3)-x(1))-min(Q_anode*x(3)-x(1)));
Results.Q_cathode_SP = (min(Q_fullcell)-min(Q_cathode*x(4)-x(2)))/(max(Q_cathode*x(4)-x(2))-min(Q_cathode*x(4)-x(2)));
Results.Q_cathode_EP = (max(Q_fullcell)-min(Q_cathode*x(4)-x(2)))/(max(Q_cathode*x(4)-x(2))-min(Q_cathode*x(4)-x(2)));

% Goodness of Fit 
Results.NRMSE = goodnessOfFit(Results.OCV_fullcell_fit',Results.OCV_fullcell','NRMSE');
Results.r2 = 1 - (norm(Results.OCV_fullcell_fit-Results.OCV_fullcell)/norm(Results.OCV_fullcell-mean(Results.OCV_fullcell))).^2;

% Recalculate fullcell OCV from balanced and aligned electrode OCVs
u_fullcell_fit = interp1(Q_cathode*x(4)-x(2), OCV_cathode, Q_fullcell,'linear')-interp1(Q_anode*x(3)-x(1), OCV_anode, Q_fullcell,'linear'); % Calculate fullcell OCV
u_dot_fullcell_fit = diff(u_fullcell_fit)./diff(Q_fullcell); % Calculate first derivative of fullcell OCV
u_rel_dot_fullcell_fit = diff(u_fullcell_fit)./diff(Q_fullcell/max(Q_fullcell)); % Calculate normalized first derivative of fullcell OCV
q_rel_dot_fullcell_fit = (Q_fullcell(2:length(Q_fullcell))/max(Q_fullcell)); % Calculate state-of-charge of fullcell OCV

% Recalculate halfcell OCVs from balanced and aligned electrode OCVs
OCV_anode_fit = interp1(Q_anode*x(3)-x(1), OCV_anode, Q_fullcell,'linear'); % Shift and scale anode OCV
u_dot_anode_fit = diff(OCV_anode_fit)./diff(Q_fullcell); % Calculate first derivative of anode OCV
u_rel_dot_anode_fit = diff(OCV_anode_fit)./diff(Q_fullcell/max(Q_fullcell)); % Calculate state-of-charge of anode OCV
OCV_cathode_fit = interp1(Q_cathode*x(4)-x(2), OCV_cathode, Q_fullcell,'linear'); % Shift and scale cathode OCV
u_dot_cathode_fit = diff(OCV_cathode_fit)./diff(Q_fullcell); % Calculate first derivative of cathode OCV
u_rel_dot_cathode_fit = diff(OCV_cathode_fit)./diff(Q_fullcell/max(Q_fullcell)); % Calculate state-of-charge of cathode OCV

% Fitting results
Results.q_dot_fullcell = q_dot_fullcell;
Results.u_dot_fullcell = u_dot_fullcell;
Results.u_dot_fullcell_fit = u_dot_fullcell_fit;
Results.u_dot_anode_fit = u_dot_anode_fit;
Results.u_dot_cathode_fit = u_dot_cathode_fit;

Results.q_rel_dot_fullcell = q_rel_dot_fullcell;
Results.u_rel_dot_fullcell = u_rel_dot_fullcell;
Results.u_rel_dot_fullcell_fit = u_rel_dot_fullcell_fit;
Results.u_rel_dot_anode_fit = u_rel_dot_anode_fit;
Results.u_rel_dot_cathode_fit = u_rel_dot_cathode_fit;

% Calculate OCV goodness of fit for 0-100% SOC range
OCV_NRMSE = goodnessOfFit(u_fullcell_fit',OCV_fullcell','NRMSE'); % Calculate OCV RMSE
OCV_r2 = 1 - (norm(u_fullcell_fit-OCV_fullcell)/norm(OCV_fullcell-mean(OCV_fullcell))).^2; % Calculate OCV R square

% Calculate DVA goodness of fit for 10-80% fast charging range
indbeg = find(q_dot_fullcell>DVAstart*max(q_dot_fullcell), 1,'first'); % 10% SOC index
indend = find(q_dot_fullcell>DVAend*max(q_dot_fullcell), 1,'first'); % 80% SOC index
DVA_NRMSE = goodnessOfFit(u_dot_fullcell_fit(indbeg:indend)',u_dot_fullcell(indbeg:indend)','NRMSE'); % Calculate DVA RMSE
DVA_r2 = 1 - (norm(u_dot_fullcell_fit(indbeg:indend)-u_dot_fullcell(indbeg:indend))/norm(u_dot_fullcell(indbeg:indend)-mean(u_dot_fullcell(indbeg:indend)))).^2; % Calculate DVA R square

% Visualize data
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color','w');

% Plot 1: OCV over charge throughput
subplot(2,1,1);
hold on;
title('Halfcell GITT aligned and balanced to fullcell GITT - DVA objective function')
xlabel('Q in As') 
ylabel('U in V')
plot(Results.Q_fullcell,Results.OCV_fullcell,'k');
plot(Results.Q_fullcell,Results.OCV_fullcell_fit,'g');
plot(Results.Q_anode,Results.OCV_anode,'k--');
plot(Results.Q_cathode,Results.OCV_cathode,'k--');
plot(Results.Q_fullcell,Results.OCV_anode_fit,'b');
plot(Results.Q_fullcell,Results.OCV_cathode_fit,'r');
legend({'Fullcell','Fullcell fit','Unused anode capacity','Unused cathode capacity','Anode','Cathode'},'interpreter','latex','Location','East');
xlim([-1500 11500]);
NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(NE(1),NE(2),{['NRMSE = ' num2str(OCV_NRMSE,4)], ['R2 = ' num2str(OCV_r2,4)]}, 'VerticalAlignment','top', 'HorizontalAlignment','right');

% Plot 2: DVA over charge throughput
subplot(2,1,2);
plot(q_dot_fullcell, u_dot_fullcell, 'k', q_dot_fullcell, u_dot_fullcell_fit,'g', q_dot_fullcell, u_dot_anode_fit, 'b', q_dot_fullcell, u_dot_cathode_fit, 'r');
hold on; 
ylabel('dU/dQ in V/As');
xlabel('Q in As');
title('DVA');
ylim([-2.5e-4 3.5e-4]);
xlim([-1500 11500]);
legend({'Fullcell','Fullcell Fit','Anode','Cathode'},'interpreter','latex','Location','East');
NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
text(NE(1),NE(2),{['NRMSE = ' num2str(DVA_NRMSE,4)], ['R2 = ' num2str(DVA_r2,4)]}, 'VerticalAlignment','top', 'HorizontalAlignment','right');

% Save results
clearvars -except Fullcell Anode Cathode Results DVA_NRMSE DVA_r2 OCV_NRMSE OCV_r2 C100_DVA_r2
save B_A Results DVA_NRMSE DVA_r2 OCV_NRMSE OCV_r2