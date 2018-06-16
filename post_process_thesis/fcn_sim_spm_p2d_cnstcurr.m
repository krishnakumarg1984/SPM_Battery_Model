function fcn_sim_spm_p2d_cnstcurr(idx,C_rate,ha,line_colors)

run('custom_colors.m');
%% Simulate cnst current SPM followed by p2d
run('user_inputs_for_sim_cnst_curr');
run('disc_spm_core_sim_cnst_curr');   % spm simulation loop is in this script
run('single_shot_p2d_cnst.m'); % p2d simulation loop is in this script

t_end_max = max(spm_sim_time_vector(end),time_vector_p2d(end));
t_end_max_order = floor(log10(t_end_max));
t_end_max_plot = ceil((t_end_max/10^t_end_max_order)/0.5)*0.5*10^t_end_max_order;

t_end_common = min(spm_sim_time_vector(end),time_vector_p2d(end));
t_common_vector = (0:Ts:t_end_common)';
max_idx = length(t_common_vector);

annotation_x = 0.025;
annotation_y = 0.225;

soc_error_pct = soc_pct_results_p2d(1:max_idx) - soc_pct_results_spm(1:max_idx);

v_error_vector = cell_voltage_results_p2d(1:max_idx) - v_cell_sim_results_spm(1:max_idx);
v_error_pct = v_error_vector*100./cell_voltage_results_p2d(1:max_idx);

abs_max_error_voltage = max(abs(v_error_pct))
% mean_error_voltage = mean(v_error_pct)
mae_voltage = mean(abs(v_error_pct))

abs_max_error_soc = max(abs(soc_error_pct))
% mean_error_soc = mean(soc_error_pct)
mae_soc = mean(abs(soc_error_pct))

% %% plot
% axes(ha(2*idx-1));
% 
% plot(time_vector_p2d,soc_pct_results_p2d,'color',line_colors(2,:),'linewidth',4);
% % plot(time_vector_p2d,cell_voltage_results_p2d,'color',line_colors(2,:)); 
% hold on;
% plot(spm_sim_time_vector,soc_pct_results_spm,'color',line_colors(1,:),'linestyle','--');
% % plot(spm_sim_time_vector,v_cell_sim_results_spm,'color',line_colors(1,:),'linestyle','--');
% hold off;
% ylim([0 100]); % for soc
% % ylim([spm_params.CutoffVoltage spm_params.CutoverVoltage]);
% ylabel('SOC (\%)');
% % ylabel('$V_\mathrm{cell}$ (V)');
% lgd = legend('P2d','SPM','location','northeast');
% lgd.Units = 'centimeters';
% ax_handle = gca;
% legend boxoff;
% ytickformat('%.0f'); % for soc
% % ytickformat('%.2f'); % for voltage
% xlim([0 t_end_max_plot]);
% no_of_x_tick_points = 5;
% no_of_y_tick_points = 5; % not required for voltage
% ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
% ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
% curtick = get(gca, 'XTick');
% set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
% % xlim([0 20000]);
% % curytick = get(gca, 'YTick');
% % set(gca, 'YTickLabel', cellstr(num2str(curytick(:)))); % remove scientific multipliers in x-axis format
% 
% % text(0.65,0.65,[num2str(C_rate(1)', '%5.2f') ' C'],'Units','normalized'); % for soc
% text(annotation_x,annotation_y,[num2str(C_rate(1)', '%5.2f') ' C'],'Units','normalized','FontWeight','bold','FontSize',13,'color',color_imp_blue); % for voltage
% 
% axes(ha(2*idx));
% plot(t_common_vector,soc_error_pct,'color',line_colors(1,:));
% % plot(t_common_vector,v_error_pct,'color',line_colors(1,:));
% ylabel('$\varepsilon_\mathrm{soc} (\%)$');
% % ylabel('$\hat{\varepsilon}_\mathrm{v} (\%)$');
% ax_handle2 = gca;
% 
% ylim([0 0.35]); % for soc error
% % ylim([-20 0]); % for voltage error
% xlim([0 t_end_max_plot]);
% ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
% ax_handle2.XAxis.TickValues = linspace(ax_handle2.XAxis.Limits(1),ax_handle2.XAxis.Limits(2),no_of_x_tick_points);
% curtick = get(gca, 'XTick');
% set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
% 
% text(annotation_x,annotation_y,[num2str(C_rate(1)', '%5.2f') ' C'],'Units','normalized','FontWeight','bold','FontSize',13,'color',color_imp_blue); % for voltage


return;