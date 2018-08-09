% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;
golden_ratio = 1.61803398875;
fig_width_factor = 1.25; % scaling factor wrt to text width
fig_ht_factor = 0.65;
figH_cm = 22.27184*fig_ht_factor;

set(0,'defaultlinelinewidth',1.5,'defaultpatchlinewidth',0.5,'defaultaxeslinewidth',1);
fig_h = clf;
fig_h.Units = 'centimeters';
figW_cm = figH_cm*fig_width_factor/golden_ratio;
% figH_cm = figW_cm/golden_ratio;
fig_h.Position = [fig_h.Position(1),fig_h.Position(2),figW_cm,figH_cm];
movegui('center');

axis_w_factor = 0.75;
axis_h_factor = axis_w_factor/(golden_ratio*1.2);


% run('setup_line_colors.m'); % deletes axes/clears plots
cbrewerintergray = [189,189,189]/255;
cbrewerdarkgray = [99,99,99]/255;
cbrewer_Gnbubu_blue = [67,162,202]/255;
% return;
%% Post-process and plot
no_of_x_tick_points = 5;
no_of_y_tick_points = 5;

clc;

load('p2d_sim_Aug_07_2018_14_50_00.mat'); % udds with 30 nodes (p2d)
load('cts_sim_Aug_08_2018_16_40_43.mat'); % udds (spm);
load('phie_sysid_Aug_07_2018_16_13_28.mat'); % udds electrolyte overpotential
pos1 = [0.2 0.575 axis_w_factor axis_h_factor]; % left, bottom, width, height
subplot('Position',pos1);
plot(time_vector_p2d,cell_voltage_results_p2d,'color',cbrewerdarkgray);
hold on;
max_valid_idx_with_phie = min([length(v_cell_sim_results_spm),length(phie_op_vector_results)]);
phie_op_vector_results = [phie_op_vector_results(1) phie_op_vector_results];
v_cell_spm_with_phie = v_cell_sim_results_spm + 0.625*(phie_op_vector_results(1:max_valid_idx_with_phie))';
plot(spm_sim_time_vector,v_cell_sim_results_spm,'color','k');
plot(spm_sim_time_vector(1:max_valid_idx_with_phie),v_cell_spm_with_phie,'color',cbrewer_Gnbubu_blue);
% title('1C discharge');
hold off;
ylabel({'$V_\mathrm{cell}$';'$\mathrm{(V)}$'});
lgd = legend('P2d','Basic SPM','Composite SPM','location','southeast');
legend boxoff;
% return;

t_end_common = min(spm_sim_time_vector(end),time_vector_p2d(end));
t_common = 0:Ts:t_end_common;
t_end_max = max(spm_sim_time_vector(end),time_vector_p2d(end));
t_end_max_order = floor(log10(t_end_max));
t_end_max_plot = ceil((t_end_max/10^t_end_max_order)/0.5)*0.5*10^t_end_max_order;
xlim([0 t_end_max_plot]);
ylim([3.6 3.95]);
ax_handle = gca;
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
ylh = get(gca,'ylabel');
gyl = get(ylh);                                                         % Object Information
ylp = get(ylh, 'Position');
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
% return;
%
max_valid_index = min([length(cell_voltage_results_p2d),length(v_cell_sim_results_spm)]);
basic_spm_abs_pct_error = 100*abs((cell_voltage_results_p2d(1:max_valid_index) - v_cell_sim_results_spm(1:max_valid_index))./cell_voltage_results_p2d(1:max_valid_index));
max_valid_index = min([length(cell_voltage_results_p2d),length(v_cell_spm_with_phie)]);
composite_spm_abs_pct_error = 100*abs((cell_voltage_results_p2d(1:max_valid_index) - v_cell_spm_with_phie(1:max_valid_index))./cell_voltage_results_p2d(1:max_valid_index));

pos2 = [0.2 0.0875 axis_w_factor axis_h_factor]; % left, bottom, width, height
subplot('Position',pos2);
plot(t_common,basic_spm_abs_pct_error,'color','k');
hold on;
plot(t_common,composite_spm_abs_pct_error,'color',cbrewer_Gnbubu_blue);
ylabel({'$|\hat{\varepsilon}_v|$';'$(\%)$'});

% title('1C constant current discharge');
xlim([0 t_end_max_plot]);
% ylim([0 11]);
ax_handle = gca;
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
xlabel('time (s)');
lgd = legend('Basic SPM','Composite SPM','location','northeast');
legend boxoff;
ylh = get(gca,'ylabel');
gyl = get(ylh);                                                         % Object Information
ylp = get(ylh, 'Position');
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')

return;

export_fig composite_spm_vcell_udds.pdf -q101
return;

extra_axis_options = 'xticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=0,/pgf/number format/fixed,/pgf/number format/fixed zerofill,},yticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=2, /pgf/number format/fixed, }, ylabel absolute, ylabel style={rotate=-90}';
custom_m2t_fcn('composite_spm_vcell_1C',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
close;

return;

