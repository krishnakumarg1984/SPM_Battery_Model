% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;

%% Simulate SPM followed by p2d
run('user_inputs_for_sim.m');
run('pre_process_script.m');
run('disc_spm_core_sim'); % spm simulation loop is in this script

%% Plot results and visualise comparison
close all;
figure(1);
h1 = subplot(211);
plot(spm_sim_time_vector,load_current_vector,'-');
ylabel('Current');
ylim([min(load_current_vector)-5 max(load_current_vector)+5]);

h2 = subplot(212);
plot(spm_sim_time_vector,v_cell_sim_results_spm); hold on;
plot(time_vector_p2d,cell_voltage_results_p2d); hold off;
ylim([spm_params.CutoffVoltage spm_params.CutoverVoltage]);
ylabel('Cell Voltage [V]');
legend('SPM','P2d','location','best');

linkaxes([h1 h2],'x');
xlim([spm_sim_time_vector(1) spm_sim_time_vector(end)]);
xlabel('Time [sec]');

figure(2);
h1 = subplot(211);
plot(spm_sim_time_vector,load_current_vector,'-'); 
ylabel('Current');
ylim([min(load_current_vector)-5 max(load_current_vector)+5]);

h2 = subplot(212);
plot(spm_sim_time_vector,soc_pct_results_spm);hold on;
plot(time_vector_p2d,soc_pct_results_p2d);hold off;
ylabel('SOC [%]');
legend('SPM','P2d','location','best');

linkaxes([h1 h2],'x');
xlim([spm_sim_time_vector(1) spm_sim_time_vector(end)]);
xlabel('Time [sec]');

clear h1 h2;

figure(1);shg;

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: