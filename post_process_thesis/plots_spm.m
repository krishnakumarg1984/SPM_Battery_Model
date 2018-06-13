%% Plot results
figure(1);
h1 = subplot(211);
plot(spm_sim_time_vector,load_current_vector,'-');
ylabel('Current');
ylim([min(load_current_vector)-5 max(load_current_vector)+5]);

h2 = subplot(212);
plot(spm_sim_time_vector,v_cell_sim_results_spm);
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
plot(spm_sim_time_vector,soc_pct_results_spm);
ylabel('SOC [%]');

linkaxes([h1 h2],'x');
xlim([spm_sim_time_vector(1) spm_sim_time_vector(end)]);
xlabel('Time [sec]');

clear h1 h2;
