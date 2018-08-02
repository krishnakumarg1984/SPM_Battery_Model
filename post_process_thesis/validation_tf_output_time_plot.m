clear; clc; close all; format short g;format compact;

set(0,'defaultaxesfontsize',12 ...
    ,'defaultaxeslinewidth',1 ...
    ,'defaultlinelinewidth',1 ...
    ,'defaultpatchlinewidth',1);

load('Qe_data_to_validate_sysid.mat');
load('Qe_data_to_validate_sysid.mat');
load('Qep_tf_4p3z_scaled.mat');
scale_factor = 1e3;
%% detrend
Qep_validate_trendobj = getTrend(Qep_data_validate);
Qep_validate_trendobj.OutputOffset = Qep_init_validate;
Qep_data_validate_detrended = detrend(Qep_data_validate,Qep_validate_trendobj);
Qep_data_validate_detrended.OutputName{1} = 'Qep_debiased';
Qep_data_validate_detrended.OutputData = Qep_data_validate_detrended.OutputData*scale_factor; % rescaling (needed?)

Qep_validate_trendobj = getTrend(Qep_data_validate);
Qep_validate_trendobj.OutputOffset = Qep_init_validate;
Qep_data_validate_detrended = detrend(Qep_data_validate,Qep_validate_trendobj);
Qep_data_validate_detrended.OutputName{1} = 'Qep_debiased';
Qep_data_validate_detrended.OutputData = Qep_data_validate_detrended.OutputData*scale_factor; % rescaling (needed?)

% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

%% 
golden_ratio = 1.61803398875;
% 
fig_width_factor = 1;
% % ax_width_factor = 0.75;
% ax_width_factor = 0.8;
% fig_ht_factor = 1;
% 
% n_axes_w = 2;   % how many horizontal/width-wise axes?
% n_axes_ht = 2;  % how many vertical axes?
% 
% if n_axes_ht <= n_axes_w
%     figW_cm = 15.74776*fig_width_factor;     % textwidth (cm) reported by LaTeX doc (with a scaling factor)
%     figH_cm = figW_cm*fig_ht_factor/golden_ratio;
% else
%     figH_cm = 22.27184*fig_ht_factor; % textheight (cm)  reported by LaTeX doc (with a scaling factor)
% %     figW_cm = figH_cm/golden_ratio;
%     figW_cm = 15.74776*fig_width_factor;
% end
% 
% ax_width = (figW_cm/n_axes_w)*ax_width_factor;
% ax_height = ax_width/golden_ratio;
% 
% % decision
% % left_margin : gap_w : right_margin
% gap_w_scale = 0.9;         % wrt to marg_w_left
% marg_w_right_scale = 0.3;   % wrt to marg_w_left
% 
% % bottom_margin : gap_ht : top_margin
% gap_ht_scale = 1.5;        % wrt to marg_ht_bottom
% marg_ht_top_scale = 1;     % wrt to marg_ht_bottom
% 
% %%
% marg_w_left = (figW_cm - n_axes_w*ax_width)/(1 + gap_w_scale + marg_w_right_scale);
% gap_w = marg_w_left*gap_w_scale;
% marg_w_right = marg_w_left*marg_w_right_scale;
% 
% marg_ht_bottom = (figH_cm - n_axes_ht*ax_height)/(1 + gap_ht_scale + marg_ht_top_scale);
% gap_ht = marg_ht_bottom*gap_ht_scale;
% marg_ht_top = marg_ht_bottom*marg_ht_top_scale;
% 
% run('setup_line_colors.m'); % deletes axes/clears plots
% [ha, pos] = tight_subplot_cm(n_axes_ht, n_axes_w, [gap_ht gap_w],[marg_ht_bottom marg_ht_top],[marg_w_left marg_w_right],figH_cm,figW_cm);
% % [ha, pos] = tight_subplot_cm(n_axes_ht, n_axes_w, gap_ht,marg_ht_bottom,marg_w_left,figH_cm,figW_cm);
% movegui('center');

%% Plots
cbrewerintergray = [189,189,189]/255;
cbrewerdarkgray = [99,99,99]/255;

no_of_x_tick_points = 8;
no_of_y_tick_points = 5;

% return;
% axes(ha(1));
fig1 = figure;
movegui('center');
fig1.Units = 'centimeters';
figW_cm = 15.74776*fig_width_factor;     % textwidth (cm) reported by LaTeX doc with a scaling factor
figH_cm = figW_cm/golden_ratio;
fig1.Position = [fig1.Position(1),fig1.Position(2),figW_cm,figH_cm];

plot(Qep_data_validate_detrended.OutputData*(1/scale_factor),'color','k');
Qep_debiased_validate_sim_scaled = (1/scale_factor)*lsim(tf(Qep_tfest_scaled),Qep_data_validate_detrended.InputData,time_p2d_validate,[],'zoh');
hold on;
plot(Qep_debiased_validate_sim_scaled,'color',cbrewerdarkgray,'linestyle','-');
xlim([0 2100]);
lgd1 = legend('P2D','SysID');  % 'location','northwest'
legend boxoff;
% return;
title('$\widetilde{Q}_{\mathrm{e,p}_\mathrm{val}}(t)$');
ax_handle = gca;
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curxtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curxtick(:)))); % remove scientific multipliers in x-axis format
% curytick = get(gca, 'YTick');
% set(gca, 'YTickLabel', cellstr(num2str(curytick(:)))); % remove scientific multipliers in y-axis format
xlabel('time (s)');
ylabel('$\mathrm{mol\, m}^{-2}$');
ylim([-6e-3 6e-3]);
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
% lgd1.Position(1) = lgd1.Position(1)+0.01;
% lgd1.Position(2) = lgd1.Position(2)+0.025;

[a,b] = max(abs(Qep_debiased_validate_sim_scaled - (Qep_data_validate_detrended.OutputData*(1/scale_factor))))
% return;
MagInset(fig1,-1, [1110, 1115,-2.47e-3 , -2.20e-3], [300, 920, -4.5e-3, -2.2e-3], {});
RMSE_validate_scaled = rms(Qep_data_validate_detrended.OutputData*(1/scale_factor) - Qep_debiased_validate_sim_scaled)
% return;

%% Export
% extra_axis_options = 'yticklabel style={/pgf/number format/1000 sep= , /pgf/number format/precision=2, /pgf/number format/fixed, }';
% extra_axis_options = 'yticklabel style = {}';
extra_axis_options = 'yticklabel style = {/pgf/number format/fixed, /pgf/number format/fixed zerofill},xticklabel style={/pgf/number format/1000 sep= },/pgfplots/tick scale binop=\times';
custom_m2t_fcn('p2d_sysid_val_qep',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
close;
% shg;
