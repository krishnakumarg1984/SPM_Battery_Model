% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;

run('user_inputs_for_sim.m');
run('pre_process_script.m');
run('core_p2d_loop_sim.m');

%% Save results to disk
save_foldername = ['p2d_results/', cellIdentifier, '/', load_profile_name];
if exist(save_foldername,'dir')==0
    mkdir(save_foldername);
end

% Replace decimal point chars in soc% string with 'p' (stands for point)
soc_init_pct_savestr = strrep(num2str(soc_init_pct),'.','p');

% clear A_disc B_disc; % potentially useful varibles. comment out for debugging
clear soc_init_pct C_rate_profile I_1C k num_iterations; % redundant info
clear t_local_finish t_local_start;
save([save_foldername, '/p2d_sim_', ...
      datestr(now, 'mmm_dd_yyyy_HH_MM_SS')]); % save workspace to file

%% Plot results
close all;
figure(1);
h1 = subplot(211);
plot(time_vector_p2d,load_curr_vector_p2d,'-');
ylabel('Current');
ylim([min(load_curr_vector_p2d)-5 max(load_curr_vector_p2d)+5]);

h2 = subplot(212);
plot(time_vector_p2d,cell_voltage_results_p2d,'m');
ylim([spm_params.CutoffVoltage spm_params.CutoverVoltage]);
ylabel('Cell Voltage [V]');

linkaxes([h1 h2],'x');
xlim([time_vector_p2d(1) time_vector_p2d(end)]);
xlabel('Time [sec]');


figure(2);
h1 = subplot(211);
plot(time_vector_p2d,load_curr_vector_p2d,'-');
ylabel('Current');
ylim([min(load_curr_vector_p2d)-5 max(load_curr_vector_p2d)+5]);

h2 = subplot(212);
plot(time_vector_p2d,soc_pct_results_p2d,'r');
ylabel('SOC [%]');

linkaxes([h1 h2],'x');
xlim([time_vector_p2d(1) time_vector_p2d(end)]);
xlabel('Time [sec]');

clear h1 h2;
figure(1);
shg;
