%% User-entered data
% case-sensitive string descriptive of cell to be simulated.
cellIdentifier = 'Northrop';

% string describing starting soc% and csv filename of load profile (time vs current through external circuit)
% a) 'cnst_dischg' b) 'cnst_chg' c) 'udds' etc
load_profile_name = 'cnst_dischg_soc_100_1C';
% load_profile_name = 'udds_soc_50';

% Input CSV-profile setup. Note: Offsets use a 0-base numbering system
soc_col            = 1; % The starting SOC is in this column of top row
profile_row_offset = 2; % Load profile input data begins only from this row

Ts       = 1;  % sec (how often are results needed?)
% Ts       = 5e-2; % sec (for capacity characterisation)
tf_user  = 100;  % sec (user-entered desired simulation end-time)
% Simulation might prematurely end if voltage/soc cutoffs are hit

termination_choice = 'max'; % valid choices are 'max' and 'min'
% The 'min' choice is helpful for trials. Whilst retaining the characteristics
% of the load profile, the user may do a short time trial simulation.
