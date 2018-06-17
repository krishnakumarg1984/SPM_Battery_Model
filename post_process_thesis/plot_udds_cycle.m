% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;
%% load sim results
load('udds.mat');

%% Plot setup user inputs
no_of_x_tick_points = 8;

fig_width_factor = 0.75; 
ax_width_factor = 0.8; % may need to modify suitably to prevent some of the dimensions from becoming negative

figW_cm = 15.74776*fig_width_factor; % textwidth (cm) reported by thesis doc from layouts package. Additonally include a scaling factor
golden_ratio = 1.618;
figH_cm = figW_cm/golden_ratio;
% fig_ht_factor = 1.0; % space for fig captions at bottom. Also, may need to modify suitably to prevent some of the dimensions from becoming negative
% figH_cm = 22.27184*fig_ht_factor; % textheight (cm) reported by thesis latex after leaving space for caption
% figW_cm = figH_cm/golden_ratio;

n_axes_w = 1; % how many horizontal/width-wise axes?
n_axes_ht = 1;  % how many vertical axes?

% left_margin : gap_w : right_margin
gap_w_scale = 0;         % wrt to marg_w_left
marg_w_right_scale = 0.3;   % wrt to marg_w_left

% bottom_margin : gap_ht : top_margin
gap_ht_scale = 0;        % wrt to marg_ht_bottom
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
time = 0:1:length(drivecycle_data.speed_met_per_sec)-1;
[ha, pos] = tight_subplot_cm(n_axes_ht, n_axes_w, [gap_ht gap_w],[marg_ht_bottom marg_ht_top],[marg_w_left marg_w_right],figH_cm,figW_cm);
axes(ha(1));
plot(time,drivecycle_data.speed_met_per_sec,'color',line_colors(1,:),'linewidth',1);

xlim([0 1400]);
ylabel('$\mathrm{speed}\; (\mathrm{ms}^{-1})$');
ax_handle = gca;
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
xlabel('time (s)');

%% Export to Tikz
extra_axis_options = 'xticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=0,/pgf/number format/fixed,/pgf/number format/fixed zerofill,},yticklabel style={/pgf/number format/1000 sep=,},';
custom_m2t_fcn('udds_cycle',[figW_cm,figH_cm]*10,[],false,extra_axis_options);

close;
