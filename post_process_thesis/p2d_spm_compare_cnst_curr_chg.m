% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;
golden_ratio = 1.618;

C_rate = -1;
fig_width_factor = 1;
ax_width_factor = 0.75;
fig_ht_factor = 0.75;

figW_cm = 15.74776*fig_width_factor;     % textwidth (cm) reported by LaTeX doc with a scaling factor
figH_cm = figW_cm/golden_ratio;

n_axes_w = 2; % how many horizontal/width-wise axes?
n_axes_ht = 2;  % how many vertical axes?

ax_width = (figW_cm/n_axes_w)*ax_width_factor;
ax_height = ax_width/golden_ratio;

% decision
% left_margin : gap_w : right_margin
gap_w_scale = 1.2;         % wrt to marg_w_left
marg_w_right_scale = 0.4;   % wrt to marg_w_left

% bottom_margin : gap_ht : top_margin
gap_ht_scale = 0.9375;        % wrt to marg_ht_bottom
marg_ht_top_scale = 0.3;     % wrt to marg_ht_bottom

%%
marg_w_left = (figW_cm - n_axes_w*ax_width)/(1 + gap_w_scale + marg_w_right_scale);
gap_w = marg_w_left*gap_w_scale;
marg_w_right = marg_w_left*marg_w_right_scale;

marg_ht_bottom = (figH_cm - n_axes_ht*ax_height)/(1 + gap_ht_scale + marg_ht_top_scale);
gap_ht = marg_ht_bottom*gap_ht_scale;
marg_ht_top = marg_ht_bottom*marg_ht_top_scale;

run('setup_line_colors.m'); % deletes axes/clears plots
[ha, pos] = tight_subplot_cm(n_axes_ht, n_axes_w, [gap_ht gap_w],[marg_ht_bottom marg_ht_top],[marg_w_left marg_w_right],figH_cm,figW_cm);
    
run('user_inputs_for_sim_cnst_curr_chg');
run('disc_spm_core_sim_cnst_curr');   % spm simulation loop is in this script
run('single_shot_p2d_cnst.m'); % p2d simulation loop is in this script

%% Post-process and plot
clc;
no_of_x_tick_points = 5;
no_of_y_tick_points = 5;

t_end_max = max(spm_sim_time_vector(end),time_vector_p2d(end));
t_end_max_order = floor(log10(t_end_max));
t_end_max_plot = ceil((t_end_max/10^t_end_max_order)/0.5)*0.5*10^t_end_max_order;

t_end_common = min(spm_sim_time_vector(end),time_vector_p2d(end));
t_common_vector = (0:Ts:t_end_common)';
max_idx = length(t_common_vector);

soc_error_pct = soc_pct_results_p2d(1:max_idx) - soc_pct_results_spm(1:max_idx);

v_error_vector = cell_voltage_results_p2d(1:max_idx) - v_cell_sim_results_spm(1:max_idx);
v_error_pct = v_error_vector*100./cell_voltage_results_p2d(1:max_idx);

axes(ha(1));
plot(time_vector_p2d,cell_voltage_results_p2d,'color',line_colors(2,:)); 
hold on;
plot(spm_sim_time_vector,v_cell_sim_results_spm,'color',line_colors(1,:),'linestyle','--');
hold off;
ylim([spm_params.CutoffVoltage spm_params.CutoverVoltage]);
ylabel('$V_\mathrm{cell}$ (V)');
lgd = legend('P2d','SPM','location','southeast');
lgd.Units = 'centimeters';
ax_handle = gca;
legend boxoff;
ytickformat('%.2f'); % for voltage
xlim([0 t_end_max_plot]);
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format

axes(ha(2));
plot(time_vector_p2d,soc_pct_results_p2d,'color',line_colors(2,:),'linewidth',4);
hold on;
plot(spm_sim_time_vector,soc_pct_results_spm,'color',line_colors(1,:),'linestyle','--');
hold off;
ylim([0 100]); % for soc
ylabel('SOC (\%)');
lgd = legend('P2d','SPM','location','southeast');
lgd.Units = 'centimeters';
ax_handle = gca;
legend boxoff;
ytickformat('%.0f'); % for soc
xlim([0 t_end_max_plot]);
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format

annotation_x = 0.50;
annotation_y = 0.55;

axes(ha(3));
abs_max_error_voltage_pct = max(abs(v_error_pct));
mae_voltage_pct = mean(abs(v_error_pct));
plot(t_common_vector,v_error_pct,'color',line_colors(1,:));
ylabel('$\hat{\varepsilon}_\mathrm{v} (\%)$');
ax_handle2 = gca;
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle2.XAxis.TickValues = linspace(ax_handle2.XAxis.Limits(1),ax_handle2.XAxis.Limits(2),no_of_x_tick_points);
text(annotation_x,annotation_y,['$\%\, |\hat{\varepsilon}_v|_\mathrm{max}:\, $' num2str(abs_max_error_voltage_pct, '%4.2f')],'Units','normalized'); % for voltage
text(annotation_x,annotation_y-0.15,['$\%\, \hat{\varepsilon}_{v_\mathrm{mae}}\,\,\,\,\,:\, $' num2str(mae_voltage_pct, '%4.2f')],'Units','normalized'); % for voltage

xlabel('time (s)');

axes(ha(4));
abs_max_error_soc = max(abs(soc_error_pct));
mae_soc = mean(abs(soc_error_pct));
plot(t_common_vector,soc_error_pct,'color',line_colors(1,:));
ylabel('$\varepsilon_\mathrm{soc} (\%)$');
ax_handle2 = gca;
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle2.XAxis.TickValues = linspace(ax_handle2.XAxis.Limits(1),ax_handle2.XAxis.Limits(2),no_of_x_tick_points);
text(annotation_x,annotation_y,['$\%\, |\hat{\varepsilon}_\mathrm{soc}|_\mathrm{max}:\, $' num2str(abs_max_error_soc, '%4.2f')],'Units','normalized'); % for voltage
text(annotation_x,annotation_y-0.15,['$\%\, \hat{\varepsilon}_{\mathrm{soc}_\mathrm{mae}}\,\,\,\,:\, $' num2str(mae_soc, '%4.2f')],'Units','normalized'); % for voltage

xlabel('time (s)');

%% Export to Tikz
extra_axis_options = 'xticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=0,/pgf/number format/fixed,/pgf/number format/fixed zerofill,},yticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=2, /pgf/number format/fixed, }, ylabel absolute,';
custom_m2t_fcn('const_curr_chg',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
close;
% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: