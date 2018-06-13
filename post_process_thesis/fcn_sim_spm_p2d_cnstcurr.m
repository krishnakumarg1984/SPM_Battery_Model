function fcn_sim_spm_p2d_cnstcurr(idx,C_rate,ha,line_colors)

%% Simulate cnst current SPM followed by p2d
run('user_inputs_for_sim_cnst_curr');
run('disc_spm_core_sim_cnst_curr');   % spm simulation loop is in this script
run('single_shot_p2d_cnst.m'); % p2d simulation loop is in this script

t_end_common = min(spm_sim_time_vector(end),time_vector_p2d(end));
t_common_vector = (0:Ts:t_end_common)';
max_idx = length(t_common_vector);
v_error_vector = cell_voltage_results_p2d(1:max_idx) - v_cell_sim_results_spm(1:max_idx);
v_error_pct = v_error_vector*100./cell_voltage_results_p2d(1:max_idx);

%% plot
axes(ha(2*idx-1));
plot(time_vector_p2d,cell_voltage_results_p2d,'color',line_colors(2,:)); 
hold on;
plot(spm_sim_time_vector,v_cell_sim_results_spm,'color',line_colors(1,:)); 
hold off;
ylim([spm_params.CutoffVoltage spm_params.CutoverVoltage]);
ylabel('$V_\mathrm{cell}$ (V)');
lgd = legend('P2d','SPM','location','northeast');
lgd.Units = 'centimeters';
ax_handle = gca;
lgd.Position(1) = ax_handle.Position(1) + 2;
lgd.Position(2) = ax_handle.Position(2) + 0.1;
legend boxoff;
ytickformat('%.2f');
no_of_x_tick_points = 5;
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
% xlim([0 20000]);
text(0.375,0.45,[num2str(C_rate(1)', '%5.2f') ' C'],'Units','normalized');


axes(ha(2*idx));
plot(t_common_vector,v_error_pct,'color',line_colors(1,:),'linestyle','--');
ylabel('$\varepsilon_\mathrm{v} (\%)$');
ax_handle2 = gca;

ylim([-20 0]);
ax_handle2.XAxis.TickValues = linspace(ax_handle2.XAxis.Limits(1),ax_handle2.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format

text(0.375,0.45,[num2str(C_rate(1)', '%5.2f') ' C'],'Units','normalized');


return;