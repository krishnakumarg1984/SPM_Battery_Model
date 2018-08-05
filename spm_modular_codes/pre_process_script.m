%% Pre-Process user data
profile_filename  = [load_profile_name,'.csv'];

% Note: a positive C-rate implies discharge and vice-versa for charge
try
    C_rate_profile = csvread(profile_filename,profile_row_offset,0);
catch
    error('a) Error in specified file, OR  b)the load profile is not in PATH. Quitting simulation ...');
end

% Compute expected end-time for allocation of storage & maximum loop indices
if strcmp(termination_choice,'max')
    t_finish = max(tf_user,C_rate_profile(end,1)); % longer of the two prevails
    % If the last time-entry in the input csv file is shorter than user-entered
    % value, then the last C-rate from the csv file is held for rest of the
    % simulation.
elseif strcmp(termination_choice,'min')
    t_finish = min(tf_user,C_rate_profile(end,1)); % longer of the two prevails
else
    error("Invalid termination choice. Valid strings are: 'max' or 'min'.");
end

% Starting SoC percentage
soc_init_pct = csvread(profile_filename,0,soc_col,[0 soc_col 0 soc_col]);

% struct of cell parameters
spm_params = parameters_spm_basic(soc_init_pct,cellIdentifier);

I_1C = spm_params.I_1C;
clear tf_user profile_row_offset soc_col profile_filename termination_choice;