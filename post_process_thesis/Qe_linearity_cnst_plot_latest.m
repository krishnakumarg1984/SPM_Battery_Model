% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;
golden_ratio = 1.61803398875;

fig_width_factor = 1;
ax_width_factor = 0.8;
fig_ht_factor = 0.65;

n_axes_w = 2;   % how many horizontal/width-wise axes?
n_axes_ht = 3;  % how many vertical axes?

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
gap_w_scale = 1.3;         % wrt to marg_w_left
marg_w_right_scale = 0.3;   % wrt to marg_w_left

% bottom_margin : gap_ht : top_margin
gap_ht_scale = 1.5;        % wrt to marg_ht_bottom
marg_ht_top_scale = 0.675;     % wrt to marg_ht_bottom

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
[ha, pos] = tight_subplot_cm(n_axes_ht, n_axes_w, [gap_ht gap_w],[marg_ht_bottom marg_ht_top],[marg_w_left marg_w_right],figH_cm,figW_cm);
% [ha, pos] = tight_subplot_cm(n_axes_ht, n_axes_w, gap_ht,marg_ht_bottom,marg_w_left,figH_cm,figW_cm);
movegui('center');
% return;
%% Post-process and plot
no_of_x_tick_points = 6;
no_of_y_tick_points = 5;

clc;
% close all;
load('results_id_Qe_vs_C_rates.mat');
C_rate1 = 1;
C_rate2 = 0.2;

alpha = 1;
beta = -2;
set(0,'defaultaxesfontsize',9,'defaultaxeslinewidth',0.5,'defaultlinelinewidth',2,'defaultpatchlinewidth',2);

%% test linearity relationship
[~, first_idx] = min(abs(C_rates_vector-C_rate1));
u1 = load_currents_vector(first_idx);
t1_vec = time_matrix(first_idx,:)';
t1_vec_valid = t1_vec(~isnan(t1_vec));
Qep1 = Qep_matrix(first_idx,1:length(t1_vec_valid))';
Qen1 = Qen_matrix(first_idx,1:length(t1_vec_valid))';

[~, second_idx] = min(abs(C_rates_vector-C_rate2));
u2 = load_currents_vector(second_idx);
t2_vec = time_matrix(second_idx,:)';
t2_vec_valid = t2_vec(~isnan(t2_vec));
Qep2 = Qep_matrix(second_idx,1:length(t2_vec_valid))';
Qen2 = Qen_matrix(second_idx,1:length(t2_vec_valid))';

C_rate3 = alpha*C_rate1 + beta*C_rate2;
[~, third_idx] = min(abs(C_rates_vector-C_rate3));
u3 = load_currents_vector(third_idx);
t3_vec = time_matrix(third_idx,:)';
t3_vec_valid = t3_vec(~isnan(t3_vec));
Qep3 = Qep_matrix(third_idx,1:length(t3_vec_valid))';
Qen3 = Qen_matrix(third_idx,1:length(t3_vec_valid))';

t_end = min([t1_vec_valid(end),t2_vec_valid(end),t3_vec_valid(end)]);

t_vec = (t0:Ts:t_end)';
last_common_index = length(t_vec);

Qep1_common = Qep1(1:last_common_index) - Qep_init;
Qep2_common = Qep2(1:last_common_index) - Qep_init;
Qep3_common = Qep3(1:last_common_index) - Qep_init;

Qen1_common = Qen1(1:last_common_index) - Qen_init;
Qen2_common = Qen2(1:last_common_index) - Qen_init;
Qen3_common = Qen3(1:last_common_index) - Qen_init;

clc;
% return;
%% Plots
xlim_max = 1500;
% return;
axes(ha(1));
plot(t_vec,Qen1_common,'color',line_colors(1,:)); 
hold on;
plot(t_vec,Qen2_common,'color',cbrewerintergray);
hold off;
title('$\widetilde{Q}_{\mathrm{e,n}_{1,2}}(t)$');
xlim([0 xlim_max]);
ax_handle = gca;
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curxtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curxtick(:)))); % remove scientific multipliers in x-axis format
text_x = 0.725;
text_yn1 = 0.8;
text_yn2 = 0.225;
text(text_x,text_yn1,'$u_1 = 60A$','Units','normalized');
text(text_x,text_yn2,'$u_2 = 12A$','Units','normalized');
% legend('$u_1 = 60A$','$u_2 = 12A$','location','northwest');
% legend boxoff;
% return;

axes(ha(2));
plot(t_vec,Qep1_common,'color',line_colors(1,:)); 
hold on;
plot(t_vec,Qep2_common,'color',cbrewerintergray);
hold off;
title('$\widetilde{Q}_{\mathrm{e,p}_{1,2}}(t)$');
xlim([0 xlim_max]);
ax_handle = gca;
ylim([ax_handle.YAxis.Limits(1) 0]);
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curxtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curxtick(:)))); % remove scientific multipliers in x-axis format
text_yp1 = 0.16;
text_yp2 = 0.75;
text(text_x,text_yp1,'$u_1 = 60A$','Units','normalized');
text(text_x,text_yp2,'$u_2 = 12A$','Units','normalized');
% return;

axes(ha(3));
plot(t_vec,Qen3_common,'color',cbrewerdarkgray);
title('$\widetilde{Q}_{\mathrm{e,n}_3}(t)$');
xlim([0 xlim_max]);
ax_handle = gca;
ylim([0 ax_handle.YAxis.Limits(2)]);
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curxtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curxtick(:)))); % remove scientific multipliers in x-axis format
% curytick = get(gca, 'YTick');
% set(gca, 'YTickLabel', cellstr(num2str(curytick(:)))); % remove scientific multipliers in x-axis format
t_Qen3 = text(text_x,text_yn1-0.075,'$u_3 = 36A$','Units','normalized');
t_Qen3.Rotation = 5;
% return;

axes(ha(4));
plot(t_vec,Qep3_common,'color',cbrewerdarkgray);
title('$\widetilde{Q}_{\mathrm{e,p}_3}(t)$');
xlim([0 xlim_max]);
ax_handle = gca;
ylim([ax_handle.YAxis.Limits(1) 0]);
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curxtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curxtick(:)))); % remove scientific multipliers in x-axis format
% curytick = get(gca, 'YTick');
% set(gca, 'YTickLabel', cellstr(num2str(curytick(:)))); % remove scientific multipliers in x-axis format
t_Qep3 = text(text_x,text_yn1-0.4,'$u_3 = 36A$','Units','normalized');
t_Qep3.Rotation = -5;
% return;

axes(ha(5));
plot(t_vec,Qen3_common,'color',line_colors(1,:));
hold on;
plot(t_vec,alpha*Qen1_common + beta*Qen2_common,'color',line_colors(1,:),'linestyle','--');
xlim([0 xlim_max]);
% title('$\qquad \qquad \widetilde{Q}_{\mathrm{e,n}_3} - \alpha \widetilde{Q}_{\mathrm{e,n}_1} - \beta \widetilde{Q}_{\mathrm{e,n}_2}$','FontSize',8);
ax_handle = gca;
ylim([0 ax_handle.YAxis.Limits(2)]);
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
xlabel('time (s)');
h_leg = legend({'$\widetilde{Q}_{\mathrm{e,n}_3}$','$\alpha\, \widetilde{Q}_{\mathrm{e,n}_1} + \beta\, \widetilde{Q}_{\mathrm{e,n}_2}$'},'location','southeast','FontSize',10);
legend boxoff;
HeightScaleFactor = 1.2;
NewHeight = h_leg.Position(4) * HeightScaleFactor;
h_leg.Position(2) = h_leg.Position(2) - (NewHeight - h_leg.Position(4));
h_leg.Position(4) = NewHeight;
% return;

axes(ha(6));
plot(t_vec,Qep3_common,'color',line_colors(1,:));
hold on;
plot(t_vec,alpha*Qep1_common + beta*Qep2_common,'color',line_colors(1,:),'linestyle','--');
xlim([0 xlim_max]);
ax_handle = gca;
ylim([ax_handle.YAxis.Limits(1) 0]);
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
xlabel('time (s)');
h_leg = legend({'$\widetilde{Q}_{\mathrm{e,p}_3}$','$\alpha\, \widetilde{Q}_{\mathrm{e,p}_1} + \beta\, \widetilde{Q}_{\mathrm{e,p}_2}$'},'location','northeast','FontSize',10);
legend boxoff;
HeightScaleFactor = 1.2;
NewHeight = h_leg.Position(4) * HeightScaleFactor;
h_leg.Position(2) = h_leg.Position(2) - (NewHeight - h_leg.Position(4));
h_leg.Position(4) = NewHeight;

return;

export_fig linearity_proof.pdf -q101
%%
return;
extra_axis_options = 'xticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=0,/pgf/number format/fixed,/pgf/number format/fixed zerofill,},yticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=2, /pgf/number format/fixed, }, ylabel absolute,';
% extra_axis_options = [];
custom_m2t_fcn('linearity_proof',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
% close;

% return;
%
