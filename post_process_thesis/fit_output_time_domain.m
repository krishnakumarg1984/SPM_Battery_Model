clear; clc; close all; format short g;format compact;

set(0,'defaultaxesfontsize',12 ...
    ,'defaultaxeslinewidth',1 ...
    ,'defaultlinelinewidth',1 ...
    ,'defaultpatchlinewidth',1);

load('Qe_data_to_fit_sysid.mat');
load('Qe_data_to_validate_sysid.mat');
load('Qen_tf_4p3z_scaled.mat');
scale_factor = 1e3;
%% detrend
Qen_fit_trendobj = getTrend(Qen_data_fit);
Qen_fit_trendobj.OutputOffset = Qen_init_fit;
Qen_data_fit_detrended = detrend(Qen_data_fit,Qen_fit_trendobj);
Qen_data_fit_detrended.OutputName{1} = 'Qen_debiased';
Qen_data_fit_detrended.OutputData = Qen_data_fit_detrended.OutputData*scale_factor; % rescaling (needed?)

Qen_validate_trendobj = getTrend(Qen_data_validate);
Qen_validate_trendobj.OutputOffset = Qen_init_validate;
Qen_data_validate_detrended = detrend(Qen_data_validate,Qen_validate_trendobj);
Qen_data_validate_detrended.OutputName{1} = 'Qen_debiased';
Qen_data_validate_detrended.OutputData = Qen_data_validate_detrended.OutputData*scale_factor; % rescaling (needed?)

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

no_of_x_tick_points = 5;
no_of_y_tick_points = 5;

% return;
% axes(ha(1));
fig1 = figure;
movegui('center');
fig1.Units = 'centimeters';
figW_cm = 15.74776*fig_width_factor;     % textwidth (cm) reported by LaTeX doc with a scaling factor
figH_cm = figW_cm/golden_ratio;
fig1.Position = [fig1.Position(1),fig1.Position(2),figW_cm,figH_cm];

plot(Qen_data_fit_detrended.OutputData*(1/scale_factor),'color','k');
Qen_debiased_fit_sim_scaled = (1/scale_factor)*lsim(tf(Qen_tfest_scaled),Qen_data_fit_detrended.InputData,time_p2d_fit,[],'zoh');
hold on;
plot(Qen_debiased_fit_sim_scaled,'color',cbrewerdarkgray,'linestyle','-');

lgd1 = legend('P2D','SysID');  % 'location','northwest'
% lgd1.Position(1) = lgd1.Position(1)-0.01;
% lgd1.Position(2) = lgd1.Position(2)+0.0375;
legend boxoff;
% return;
title('$\widetilde{Q}_{\mathrm{e,n}_\mathrm{train}}(t)$');
xlim([0 1200]);
ax_handle = gca;
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curxtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curxtick(:)))); % remove scientific multipliers in x-axis format
% curytick = get(gca, 'YTick');
% set(gca, 'YTickLabel', cellstr(num2str(curytick(:)))); % remove scientific multipliers in y-axis format
xlabel('time (s)');
ylabel('$\mathrm{mol\, m}^{-2}$');
InSet = get(gca, 'TightInset');
set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);

MagInset(fig1,-1, [288, 300,-1.175e-3 , -0.9e-3], [400, 750, -4.5e-3, -1.75e-3], {});
RMSE_fit_scaled = rms(Qen_data_fit_detrended.OutputData*(1/scale_factor) - Qen_debiased_fit_sim_scaled)
return;

%% Export
% extra_axis_options = 'yticklabel style={/pgf/number format/1000 sep= , /pgf/number format/precision=2, /pgf/number format/fixed, }';
% extra_axis_options = 'yticklabel style = {}';
extra_axis_options = 'yticklabel style = {/pgf/number format/fixed, /pgf/number format/fixed zerofill},xticklabel style={/pgf/number format/1000 sep= },/pgfplots/tick scale binop=\times';
custom_m2t_fcn('p2d_sysid_train_qen',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
close;
% shg;
