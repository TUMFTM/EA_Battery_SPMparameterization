%% Compute Initial Solid Concentrations from SOC

function SOC = soc_from_cell_state(p, soc_p, soc_n)
SOC_p_temp = ((1-soc_p) - p.B_A.Q_cathode_SP)/(p.B_A.Q_cathode_EP - p.B_A.Q_cathode_SP)*100;
SOC_n_temp = (soc_n-p.B_A.Q_anode_SP)/(p.B_A.Q_anode_EP - p.B_A.Q_anode_SP)*100;


SOC = (SOC_p_temp+SOC_n_temp)/2;

% some calculation error leads to tiny negative SOC at empty Cell
if SOC<0
    SOC=0;
end