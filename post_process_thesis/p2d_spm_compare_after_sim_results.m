% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;

fig_width_factor = 0.8; 
ax_width_factor = 0.75; % may need to modify suitably to prevent some of the dimensions from becoming negative

figW_cm = 15.74776*fig_width_factor; % textwidth (cm) reported by thesis doc from layouts package. Additonally include a scaling factor
golden_ratio = 1.618;
% figH_cm = figW_cm/golden_ratio;
fig_ht_factor = 0.85; % may need to modify suitably to prevent some of the dimensions from becoming negative
figH_cm = 22.27184*fig_ht_factor; % textheight (cm) reported by thesis latex after leaving space for caption

n_axes_w = 1; % how many horizontal/width-wise axes?
n_axes_ht = 3;  % how many vertical axes?

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

%% load sim results
load('p2d_sim_Jun_16_2018_23_47_38.mat'); % udds p2d
load('disc_sim_Jun_16_2018_23_43_22.mat'); % udds spm

t_end_max = max(spm_sim_time_vector(end),time_vector_p2d(end));
t_end_max_order = floor(log10(t_end_max));
t_end_max_plot = ceil((t_end_max/10^t_end_max_order)/0.5)*0.5*10^t_end_max_order;

t_end_common = min(spm_sim_time_vector(end),time_vector_p2d(end));
t_common_vector = (0:Ts:t_end_common)';
max_idx = length(t_common_vector);

