% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;

%% Fit data
load('Qe_data_to_fit_sysid.mat');
u_fit = Qen_data_fit.InputData;
u_fit = u_fit - mean(u_fit);
y1_fit = Qen_data_fit.OutputData;
y1_fit = y1_fit - mean(y1_fit);
% subplot(211);
% plot(u_fit);
% subplot(212);
% plot(y_fit);
Ts = 1;
z1_fit = iddata(y1_fit,u_fit,Ts);
[m1,p1,w1]=bode(spa(z1_fit));
w1=squeeze(w1); m1=squeeze(m1);

% NN = struc(1:4,1:3,0); % na, nb, nk
% ze = z1_fit(1:600);
% zv = z1_fit(601:1200);
% V = ivstruc(ze,zv,NN);
% order = selstruc(V);
% return;


load('Qe_data_to_validate_sysid.mat');
u_validate = Qen_data_validate.InputData;
u_validate = u_validate - mean(u_validate);
y2_validate = Qep_data_validate.OutputData;
y2_validate = y2_validate - mean(y2_validate);
z2_validate = iddata(y2_validate,u_validate,Ts);
[m2,p2,w2]=bode(spa(z2_validate));
w2=squeeze(w2); m2=squeeze(m2);

% return;
% u_val = Qen_data_validate.InputData;
% u_val = u_val - mean(u_val);
%% 
golden_ratio = 1.61803398875;
% 
fig_width_factor = 1;
ax_width_factor = 0.7;
fig_ht_factor = 0.7;

n_axes_w = 2;   % how many horizontal/width-wise axes?
n_axes_ht = 1;  % how many vertical axes?

if n_axes_ht <= n_axes_w
    figW_cm = 15.74776*fig_width_factor;     % textwidth (cm) reported by LaTeX doc (with a scaling factor)
    figH_cm = figW_cm*fig_ht_factor/golden_ratio;
else
    figH_cm = 22.27184*fig_ht_factor; % textheight (cm)  reported by LaTeX doc (with a scaling factor)
%     figW_cm = figH_cm/golden_ratio;
    figW_cm = 15.74776*fig_width_factor;
end

ax_width = (figW_cm/n_axes_w)*ax_width_factor;
ax_height = ax_width/golden_ratio;

% decision
% left_margin : gap_w : right_margin
gap_w_scale = 1;         % wrt to marg_w_left
marg_w_right_scale = 0.3;   % wrt to marg_w_left

% bottom_margin : gap_ht : top_margin
gap_ht_scale = 1;        % wrt to marg_ht_bottom
marg_ht_top_scale = 0.7;     % wrt to marg_ht_bottom

%%
marg_w_left = (figW_cm - n_axes_w*ax_width)/(1 + gap_w_scale + marg_w_right_scale);
gap_w = marg_w_left*gap_w_scale;
marg_w_right = marg_w_left*marg_w_right_scale;

marg_ht_bottom = (figH_cm - n_axes_ht*ax_height)/(1 + gap_ht_scale + marg_ht_top_scale);
gap_ht = marg_ht_bottom*gap_ht_scale;
marg_ht_top = marg_ht_bottom*marg_ht_top_scale;

run('setup_line_colors.m'); % deletes axes/clears plots
cbrewerintergray = [189,189,189]/255;
cbrewerdarkgray = [99,99,99]/255;
% return;

%% Post-process and plot
no_of_x_tick_points = 6;
no_of_y_tick_points = 6;

clc;

%% Plots
% fig_h = clf;
% fig_h.Units = 'centimeters';
% figW_cm = 15.74776*fig_width_factor;     % textwidth (cm) reported by LaTeX doc with a scaling factor
% figH_cm = figW_cm/golden_ratio;
set(0,'defaultaxesfontsize',12,'defaultlinelinewidth',2,'defaultpatchlinewidth',0.5,'defaultaxeslinewidth',1);
% fig_h.Position = [fig_h.Position(1),fig_h.Position(2),figW_cm,figH_cm];
[ha, pos] = tight_subplot_cm(n_axes_ht, n_axes_w, [gap_ht gap_w],[marg_ht_bottom marg_ht_top],[marg_w_left marg_w_right],figH_cm,figW_cm);
movegui('center');
axes(ha(1));
semilogx(w1,20 * log10(m1),'color',line_colors(1,:)); hold on;
load('Qen_tf_4p3z_scaled.mat');
[m_Tf,p_tf,w_tf]=bode(Qen_tfest_scaled);
w_tf=squeeze(w_tf); m_Tf=squeeze(m_Tf);
semilogx(w_tf,20 * log10(m_Tf/1000),'color',cbrewerdarkgray);
xlabel('Frequency, $\omega$ (rad/s)');
ylabel('dB');
title('$|\widetilde{Q}_{\mathrm{e,n}_\mathrm{train}}| $');
% return;

xlim([1e-2 1e1]);
xticks([1e-2 1e-1 1 10]);
ax_handle = gca;
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
% ax_handle.XAxis.TickValues = logspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);

axes(ha(2));
semilogx(w2,20 * log10(m2),'color',line_colors(1,:));hold on;
load('Qep_tf_4p3z_scaled.mat');
[m_Tf_val,p_tf_val,w_tf_val]=bode(Qep_tfest_scaled);
w_tf_val=squeeze(w_tf_val); m_Tf_val=squeeze(m_Tf_val);
semilogx(w_tf_val,20 * log10(m_Tf_val/1e4),'color',cbrewerdarkgray);

xlabel('Frequency, $\omega$ (rad/s)');
ylabel('dB');
title('$|\widetilde{Q}_{\mathrm{e,p}_\mathrm{val}}| $');
xlim([1e-2 1e1]);
xticks([1e-2 1e-1 1 10]);
ax_handle = gca;
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage

% export_fig linearity_proof.pdf -q101
%%
return;
extra_axis_options = 'ylabel absolute,';
% extra_axis_options = [];
custom_m2t_fcn('bode_mag',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
% close;

return;
