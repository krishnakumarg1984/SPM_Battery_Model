% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

% clear;clc; format short g; format compact; close all;
function p2d_spm_compare()
%% Simulate SPM followed by p2d
run('user_inputs_for_sim.m');
run('pre_process_script.m');
run('disc_spm_core_sim');   % spm simulation loop is in this script
% run('core_p2d_loop_sim.m'); % p2d simulation loop is in this script
run('single_shot_p2d_cnst.m'); % p2d simulation loop is in this script

% %% characterise error
% t_end_common = min(spm_sim_time_vector(end),time_vector_p2d(end));
% t_common = (0:Ts:t_end_common)';
% max_idx = length(t_common);
% v_error_vector = cell_voltage_results_p2d(1:max_idx) - v_cell_sim_results_spm(1:max_idx);
% max_v_error = max(abs(v_error_vector));

return;
% save('dischg_Cby100_soc100_p2d_basicspm.mat');
% %% Plot results and visualise comparison
% close all;
% figure(1);
% h1 = subplot(211);
% plot(spm_sim_time_vector,load_current_vector,'-');
% ylabel('Current, $I_\mathrm{load}$ (A)');
% ylim([min(load_current_vector)-5 max(load_current_vector)+5]);
% 
% h2 = subplot(212);
% plot(spm_sim_time_vector,v_cell_sim_results_spm); hold on;
% plot(time_vector_p2d,cell_voltage_results_p2d); hold off;
% ylim([spm_params.CutoffVoltage spm_params.CutoverVoltage]);
% ylabel('Cell Voltage, $V_\mathrm{cell}$ [V]');
% legend('SPM','P2d','location','best');
% 
% linkaxes([h1 h2],'x');
% xlim([spm_sim_time_vector(1) spm_sim_time_vector(end)]);
% xlabel('time (s)');
% 
% figure(2);
% h1 = subplot(211);
% plot(spm_sim_time_vector,load_current_vector,'-'); 
% ylabel('Current, $I_\mathrm{load}$ (A)');
% ylim([min(load_current_vector)-5 max(load_current_vector)+5]);
% 
% h2 = subplot(212);
% plot(spm_sim_time_vector,soc_pct_results_spm);hold on;
% plot(time_vector_p2d,soc_pct_results_p2d);hold off;
% ylabel('SOC (\%)');
% legend('SPM','P2d','location','best');
% 
% linkaxes([h1 h2],'x');
% xlim([spm_sim_time_vector(1) spm_sim_time_vector(end)]);
% xlabel('time (s)');
% 
% clear h1 h2;
% 
% figure(1);shg;

% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: