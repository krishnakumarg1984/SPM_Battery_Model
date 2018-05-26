% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;

%% User-entered data
% case-sensitive string descriptive of cell to be simulated.
cellIdentifier = 'Northrop';

% string describing starting soc% and csv filename of load profile (time vs current through external circuit)
load_profile_name = 'cnst_dischg'; % a) 'cnst_dischg' b) 'cnst_chg' c) 'udds' etc

% Input CSV-profile setup. Note: Offsets use a 0-base numbering system
soc_col            = 1; % The starting SOC is in this column of top row
profile_row_offset = 2; % Load profile input data begins only from this row

Ts       = 0.5; % sec (how often are results needed?)
tf_user  = 1.5; % sec (user-entered desired simulation end-time)
% Simulation might prematurely end if voltage/soc cutoffs are hit

termination_choice = 'max'; % valid choices are 'max' and 'min'
% The 'min' choice is helpful for trials. Whilst retaining the characteristics
% of the load profile, the user may do a short time trial simulation.

%% Pre-Process user data
profile_filename  = [load_profile_name,'.csv'];

% Note: a positive C-rate implies discharge and vice-versa for charge
try
    C_rate_profile = csvread(profile_filename,profile_row_offset,0);
catch
    error('Invalid load profile specified. Quitting simulation ...');
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

%% Define the State-eqn  and Output equation for simulation
stateEqn   = @spm_cts_stateEqn_three_states;
outputEqn  = @spm_three_states_battery_voltage;

%% Allocate storage for simulated quantities
num_iterations = ceil(t_finish/Ts) + 1; % max no. of steps (assuming no cutoff)

spm_sim_time_vector        = nan(num_iterations,1);
load_current_vector        = nan(num_iterations,1);
v_cell_sim_results_spm     = nan(num_iterations,1);
soc_pct_results_spm        = nan(num_iterations,1);
cs_avg_neg_sim_results_spm = nan(num_iterations,1);
q_pos_sim_results_spm      = nan(num_iterations,1);
q_neg_sim_results_spm      = nan(num_iterations,1);

%% Initialise SPM state vector and all other simulated quantities
spm_sim_time_vector(1)        = 0;
soc_pct_results_spm(1)        = soc_init_pct;
cs_avg_neg_sim_results_spm(1) = spm_params.cs_n_init;
q_pos_sim_results_spm(1)      = 0;
q_neg_sim_results_spm(1)      = 0;

% load current applied at t = t0
load_current_vector(1) = I_1C*interp1(C_rate_profile(:,1),C_rate_profile(:,2),spm_sim_time_vector(1),'previous','extrap');

x_spm_init = [q_pos_sim_results_spm(1); ...
              q_neg_sim_results_spm(1); ...
              cs_avg_neg_sim_results_spm(1)];

v_cell_sim_results_spm(1) = outputEqn(x_spm_init,load_current_vector(1),spm_params);

t_local_start      = 0;
t_local_finish     = Ts;
x_spm_local_finish = x_spm_init;

clear x_init t_finish q_pos_init q_neg_init;

%% Simulate the SPM
progressbarText(0);

for k = 2:num_iterations  % Need solution at k-th time-step
    load_current_vector(k)        = I_1C*interp1(C_rate_profile(:,1),C_rate_profile(:,2),spm_sim_time_vector(1),'previous','extrap'); % load current that was held constant from (k-1) to (k)
    t_span                        = [t_local_start t_local_finish];
    [~,x_spm_local_finish_matrix] = ode45(@(t,x_spm_local_finish) stateEqn(x_spm_local_finish,load_current_vector(k),spm_params), t_span, x_spm_local_finish);
    x_spm_local_finish            = x_spm_local_finish_matrix(end,:)';

    q_pos_sim_results_spm(k)      = x_spm_local_finish(1);
    q_neg_sim_results_spm(k)      = x_spm_local_finish(2);
    cs_avg_neg_sim_results_spm(k) = x_spm_local_finish(3);

    soc_pct_results_spm(k)        = 100*((cs_avg_neg_sim_results_spm(k)/spm_params.cs_max_n) - spm_params.theta_min_neg)/(spm_params.theta_max_neg - spm_params.theta_min_neg);
    v_cell_sim_results_spm(k)     = outputEqn(x_spm_local_finish,load_current_vector(k),spm_params);

    overall_exit_status = check_termination(soc_pct_results_spm(k),v_cell_sim_results_spm(k),spm_params);
    if overall_exit_status ~= 0 % check for violation of cut-off conditions
        k = k - 1; % Values in the last simulated index are incorrect.
        fprintf('Exiting simulation ...\n');
        break;
    end

    spm_sim_time_vector(k) = t_local_finish;
    t_local_start          = t_local_finish;
    t_local_finish         = t_local_start + Ts;
    progressbarText((k-1)/num_iterations);
end
fprintf('\n');

spm_sim_time_vector        = spm_sim_time_vector(1:k);
load_current_vector        = load_current_vector(1:k);
q_pos_sim_results_spm      = q_pos_sim_results_spm(1:k);
q_neg_sim_results_spm      = q_neg_sim_results_spm(1:k);
cs_avg_neg_sim_results_spm = cs_avg_neg_sim_results_spm(1:k);
soc_pct_results_spm        = soc_pct_results_spm(1:k);
v_cell_sim_results_spm     = v_cell_sim_results_spm(1:k);

%% Save results to disk
if exist('spm_results','dir')==0
    mkdir('spm_results');
end

% Replace decimal point chars in soc% string with 'p' (stands for point)
soc_init_pct_savestr = strrep(num2str(soc_init_pct),'.','p');

% clear A_disc B_disc; % potentially useful varibles. comment out for debugging
clear soc_init_pct C_rate_profile I_1C k num_iterations; % redundant info
clear x_spm_init x_spm_local_finish t_local_finish t_local_start;
save(['spm_results/cts_sim_', cellIdentifier, '_', load_profile_name, ...
      '_initial_soc_', soc_init_pct_savestr, 'pct', ...
      datestr(now, '_mmm_dd_yyyy_HH_MM_SS')]); % save workspace to file

%% Plot results
close all;
figure(1);
h1 = subplot(211);
plot(spm_sim_time_vector,load_current_vector,'-');
ylabel('Current');
ylim([min(load_current_vector)-5 max(load_current_vector)+5]);

h2 = subplot(212);
plot(spm_sim_time_vector,v_cell_sim_results_spm,'m');
ylim([spm_params.CutoffVoltage spm_params.CutoverVoltage]);
ylabel('Cell Voltage [V]');

linkaxes([h1 h2],'x');
xlim([spm_sim_time_vector(1) spm_sim_time_vector(end)]);
xlabel('Time [sec]');


figure(2);
h1 = subplot(211);
plot(spm_sim_time_vector,load_current_vector,'-');
ylabel('Current');
ylim([min(load_current_vector)-5 max(load_current_vector)+5]);

h2 = subplot(212);
plot(spm_sim_time_vector,soc_pct_results_spm,'r');
ylabel('SOC [%]');

linkaxes([h1 h2],'x');
xlim([spm_sim_time_vector(1) spm_sim_time_vector(end)]);
xlabel('Time [sec]');

clear h1 h2;
figure(1);
shg;

