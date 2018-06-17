% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;
%% load sim results
load('p2d_sim_Jun_16_2018_23_47_38.mat'); % udds p2d
load('disc_sim_Jun_16_2018_23_43_22.mat'); % udds spm

%% Error calcs
t_end_max = max(spm_sim_time_vector(end),time_vector_p2d(end));
t_end_max_order = floor(log10(t_end_max));
t_end_max_plot = ceil((t_end_max/10^t_end_max_order)/0.1)*0.1*10^t_end_max_order;

t_end_common = min(spm_sim_time_vector(end),time_vector_p2d(end));
t_common_vector = (0:Ts:t_end_common)';
max_idx = length(t_common_vector);

v_error_vector = cell_voltage_results_p2d(1:max_idx) - v_cell_sim_results_spm(1:max_idx);
soc_error_pct = soc_pct_results_p2d(1:max_idx) - soc_pct_results_spm(1:max_idx);

abs_max_error_voltage = max(abs(v_error_vector))*1000 %mV
abs_max_error_soc = max(abs(soc_error_pct))

mae_voltage = mean(abs(v_error_vector))*1000 %mV
mae_soc = mean(abs(soc_error_pct))

rmse_voltage = rms(v_error_vector)*1000 %mV
rmse_soc = rms(soc_error_pct)
return;
%% Plot setup user inputs
fig_width_factor = 1; 
ax_width_factor = 0.775; % may need to modify suitably to prevent some of the dimensions from becoming negative

% figW_cm = 15.74776*fig_width_factor; % textwidth (cm) reported by thesis doc from layouts package. Additonally include a scaling factor
golden_ratio = 1.618;
% figH_cm = figW_cm/golden_ratio;
fig_ht_factor = 1.0; % space for fig captions at bottom. Also, may need to modify suitably to prevent some of the dimensions from becoming negative
figH_cm = 22.27184*fig_ht_factor; % textheight (cm) reported by thesis latex after leaving space for caption
figW_cm = figH_cm/golden_ratio;

n_axes_w = 1; % how many horizontal/width-wise axes?
n_axes_ht = 3;  % how many vertical axes?

% left_margin : gap_w : right_margin
gap_w_scale = 0;         % wrt to marg_w_left
marg_w_right_scale = 0.4;   % wrt to marg_w_left

% bottom_margin : gap_ht : top_margin
gap_ht_scale = 1;        % wrt to marg_ht_bottom
marg_ht_top_scale = 0.2;     % wrt to marg_ht_bottom

%% Plot setup calcs
ax_width = (figW_cm/n_axes_w)*ax_width_factor;
ax_height = ax_width/golden_ratio;

marg_w_left = (figW_cm - n_axes_w*ax_width)/(1 + gap_w_scale + marg_w_right_scale);
gap_w = marg_w_left*gap_w_scale;
marg_w_right = marg_w_left*marg_w_right_scale;

marg_ht_bottom = (figH_cm - n_axes_ht*ax_height)/(1 + gap_ht_scale + marg_ht_top_scale);
gap_ht = marg_ht_bottom*gap_ht_scale;
marg_ht_top = marg_ht_bottom*marg_ht_top_scale;

run('setup_line_colors.m'); % deletes axes/clears plots

%% Plot
set(0,'defaultaxesfontsize',12 ...
    ,'defaultaxeslinewidth',2 ...
    ,'defaultlinelinewidth',1 ...
    ,'defaultpatchlinewidth',1);
no_of_x_tick_points = 8;
% no_of_y_tick_points = 5;
close all;clc;
[ha, fig_handle] = tight_subplot_cm(n_axes_ht, n_axes_w, [gap_ht gap_w],[marg_ht_bottom marg_ht_top],[marg_w_left marg_w_right],figH_cm,figW_cm);
% return;
fig_handle.Units = 'normalized';
curr_position = fig_handle.Position;
set(fig_handle, 'Position',[curr_position(1)+0.35   curr_position(2)+0.1   curr_position(3)   curr_position(4)]);

axes(ha(1));
plot(time_vector_p2d,load_curr_vector_p2d,'color',line_colors(1,:));
ylabel('Current (A)');
xlim([0 t_end_max_plot]);
ax_handle = gca;
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
% ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage

axes(ha(2));
plot(time_vector_p2d,cell_voltage_results_p2d,'color',line_colors(2,:));
hold on;
plot(spm_sim_time_vector,v_cell_sim_results_spm,'color',line_colors(1,:));
hold off;
ylabel('$V_\mathrm{cell}$ (V)');
lgd = legend('P2d','SPM','location','northeast');
lgd.Units = 'normalized';
% return;
% lgd.Position = [lgd.Position(1) lgd.Position(2)+0.015 lgd.Position(3) lgd.Position(4)]; 
legend boxoff;
ytickformat('%.2f');
xlim([0 t_end_max_plot]);
ax_handle = gca;
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
% ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
% MagInset(fig_handle, -1, [600,700,3.75,3.90], [850,1100,3.7,3.9], {'NW','NW';'SE','SE'});

axes(ha(3));
plot(time_vector_p2d,soc_pct_results_p2d,'color',line_colors(2,:),'linewidth',4);
hold on;
plot(spm_sim_time_vector,soc_pct_results_spm,'color',line_colors(1,:),'linestyle','--','linewidth',2);
ylabel('SOC (\%)');
lgd = legend('P2d','SPM','location','northeast');
lgd.Units = 'centimeters';
legend boxoff;
xlim([0 t_end_max_plot]);
ax_handle = gca;
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
% ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage

xlabel('time (s)');

return;
%% Export to Tikz
extra_axis_options = 'xticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=0,/pgf/number format/fixed,/pgf/number format/fixed zerofill,},yticklabel style={/pgf/number format/1000 sep=,}, ylabel absolute,';
custom_m2t_fcn('udds_I_v_soc',[figW_cm,figH_cm]*10,[],false,extra_axis_options);

close;
% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: