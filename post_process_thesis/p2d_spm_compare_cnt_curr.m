% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;
golden_ratio = 1.618;
ax_width_factor = 0.75;
fig_ht_factor = 0.75;

figW_cm = 15.74776;     % textwidth (cm) reported by LaTeX doc
figH_cm = 22.27184*fig_ht_factor; % fig height after leaving space for caption

C_rate_choices = [0.2,0.5,1,3];
n_axes_w = 2; % how many horizontal/width-wise axes?
n_axes_ht = length(C_rate_choices);  % how many vertical axes?

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

for idx = 1:length(C_rate_choices)
    C_rate = C_rate_choices(idx);
    fcn_sim_spm_p2d_cnstcurr(idx,C_rate,ha,line_colors);
end

axes(ha(2*idx-1));
xlabel('time (s)');
axes(ha(2*idx));
xlabel('time (s)');

%% Export to Tikz
extra_axis_options = 'xticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=0,/pgf/number format/fixed,/pgf/number format/fixed zerofill,},yticklabel style={/pgf/number format/1000 sep=,},';
custom_m2t_fcn('const_curr_dischg_soc',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
% custom_m2t_fcn('const_curr_dischg_voltage',[figW_cm,figH_cm]*10,[],false,extra_axis_options);

close;
% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: