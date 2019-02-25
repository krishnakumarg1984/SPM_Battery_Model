%% Define the State-eqn and Output equation for simulation
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
    load_current_vector(k)        = I_1C*interp1(C_rate_profile(:,1),C_rate_profile(:,2),spm_sim_time_vector(k-1),'previous','extrap'); % load current that was held constant from (k-1) to (k)
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
