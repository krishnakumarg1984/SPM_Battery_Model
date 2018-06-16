% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;

%% Load user data, pre-process and run simulation loop
run('user_inputs_for_sim.m');
run('pre_process_script.m');
run('disc_spm_core_sim.m');

%% Save results to disk
save_foldername = ['spm_results/', cellIdentifier, '/', load_profile_name];
if exist(save_foldername,'dir')==0
    mkdir(save_foldername);
end

% Replace decimal point chars in soc% string with 'p' (stands for point)
soc_init_pct_savestr = strrep(num2str(soc_init_pct),'.','p');

% clear A_disc B_disc; % potentially useful varibles. comment out for debugging
clear soc_init_pct C_rate_profile I_1C k num_iterations; % redundant info
clear x_spm_init x_spm_local_finish t_local_finish t_local_start;
save([save_foldername, '/disc_sim_', ...
      datestr(now, 'mmm_dd_yyyy_HH_MM_SS')]); % save workspace to file

%% Plot results
close all;
figure(1);
h1 = subplot(211);
plot(spm_sim_time_vector,load_current_vector,'-');
ylabel('Current');
ylim([min(load_current_vector)-5 max(load_current_vector)+5]);

h2 = subplot(212);
plot(spm_sim_time_vector,v_cell_sim_results_spm,'m');
% ylim([spm_params.CutoffVoltage spm_params.CutoverVoltage]);
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
ylabel('SOC [\%]');

linkaxes([h1 h2],'x');
xlim([spm_sim_time_vector(1) spm_sim_time_vector(end)]);
xlabel('Time [sec]');

clear h1 h2;
figure(1);
shg;

