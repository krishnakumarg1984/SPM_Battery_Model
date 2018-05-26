function param = parameters_spm_basic(init_cell_soc_percent, cellIdentifier)
% cellIdentifier is a string descriptive of the cell to be simulated.
% A matlab script of the filename pattern  'cell_params_cellIdentifier.m'
% needs to be present in matlab's PATH.

% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

%% Constants
param.F    = 96487;  % Faraday Constant  [C/mol]
param.R    = 8.314;  % Gas constant      [J / (mol K)]
param.Tref = 298.15; % Initial temperature of the cell [K]

%% Cut-off Conditions
param.CutoffVoltage  = 2.5; % Cutoff voltage [V]
param.CutoverVoltage = 4.3; % Cutover voltage [V]
param.CutoffSOC      = 0;   % Cutoff SOC percentage
param.CutoverSOC     = 100; % Cutover SOC percentage

%% Append cell specific parameters to the 'param' struct
try
    run(['cell_params_',cellIdentifier]); % append cell-specific params to 'params' struct
catch
    error('The parameterisation file with specified identifier does not exist in PATH');
end

%% Computation of initial concentration of Li-ions in the solid phase [mol/m^3]
param.init_cell_soc = init_cell_soc_percent/100; % convert to fraction
param.cs_p_init = ((param.init_cell_soc*(param.theta_max_pos-param.theta_min_pos) + param.theta_min_pos))*param.cs_max_p;
param.cs_n_init = ((param.init_cell_soc*(param.theta_max_neg-param.theta_min_neg) + param.theta_min_neg))*param.cs_max_n;

end
