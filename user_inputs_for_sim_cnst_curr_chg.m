%% User-entered data
% case-sensitive string descriptive of cell to be simulated.
cellIdentifier = 'Northrop';

Ts       = 0.5;          % sec (how often are results needed?)
t_finish = abs(3600/C_rate);  % sec (user-entered desired simulation end-time)
% Simulation might prematurely end if voltage/soc cutoffs are hit

soc_init_pct = 0.1; % in percentage

% struct of cell parameters
spm_params = parameters_spm_basic(soc_init_pct,cellIdentifier);
I_1C = spm_params.I_1C;