% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;
fig_width_factor = 1; % scaling factor wrt to text width
golden_ratio = 1.618;

C_rate = 1;
run('user_inputs_for_sim_cnst_curr_chg');
run('single_shot_p2d_cnst.m'); % p2d simulation loop is in this script

%% Set-up x-grid for plot
x_pos = linspace(param_p2d{1}.len_p/(2*param_p2d{1}.Np), param_p2d{1}.len_p-(param_p2d{1}.len_p/(2*param_p2d{1}.Np)), param_p2d{1}.Np);
x_sep = param_p2d{1}.len_p + linspace(param_p2d{1}.len_s/(2*param_p2d{1}.Ns), param_p2d{1}.len_s-(param_p2d{1}.len_s/(2*param_p2d{1}.Ns)), param_p2d{1}.Ns);
x_neg = param_p2d{1}.len_p + param_p2d{1}.len_s + linspace(param_p2d{1}.len_n/(2*param_p2d{1}.Nn), param_p2d{1}.len_n-(param_p2d{1}.len_n/(2*param_p2d{1}.Nn)), param_p2d{1}.Nn);
x = [x_pos, x_sep, x_neg];
x_um = x*1e6;
x_um_full_length = 1e6*linspace(0,param_p2d{1}.len_p+param_p2d{1}.len_s+param_p2d{1}.len_n,5*(param_p2d{1}.Np+param_p2d{1}.Ns+param_p2d{1}.Nn));
%% transient plots
clc;
t_plot_transient_choices = [5,10,20];
t_plot_transient_indices = find(ismember(time_vector_p2d,t_plot_transient_choices));

fig_h = clf;
fig_h.Units = 'centimeters';
figW_cm = 15.74776*fig_width_factor;     % textwidth (cm) reported by LaTeX doc with a scaling factor
figH_cm = figW_cm/golden_ratio;
set(0,'defaultlinelinewidth',1.5,'defaultpatchlinewidth',0.5,'defaultaxeslinewidth',1);
fig_h.Position = [fig_h.Position(1),fig_h.Position(2),figW_cm,figH_cm];
N = length(t_plot_transient_choices);
% OrRdmap = brewermap(N+1,'OrRd'); % OrRd upto 4 colors is a) colorblind-safe, b) print-safe and c)photocopy-safe
OrRdmap = flipud(magma(N));
OrRdmap_colors = OrRdmap; 
% OrRdmap_colors = OrRdmap(2:end,:); 
for loop_idx = 1:length(t_plot_transient_choices)
    t_plot_idx = t_plot_transient_indices(loop_idx);
    ce_p2d_t_plot_idx = ce_p2d(t_plot_idx,:);
    ce_p2d_t_plot_idx_full_length = interp1(x_um,ce_p2d_t_plot_idx,x_um_full_length,'linear','extrap');
    plot(x_um_full_length,ce_p2d_t_plot_idx_full_length,'color',OrRdmap_colors(loop_idx,:)); hold on;
end

% xlabel('$x\, \mathrm{(\mu\, m)}$');
ylabel('$c_e\, (\mathrm{mol\, m}^{-3})$');
ax_handle = gca;
InSet = get(ax_handle, 'TightInset');
set(ax_handle, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
xticks(1e6*[0 param_p2d{1}.len_p param_p2d{1}.len_p+param_p2d{1}.len_s param_p2d{1}.len_p+param_p2d{1}.len_s+param_p2d{1}.len_n]);
xticklabels({'0','$l_\mathrm{pos}$','$l_\mathrm{pos}\! + l_\mathrm{sep}$','$l_\mathrm{tot}$'})
xlim([0 x_um_full_length(end)]);


shg;

%% SS Plot
t_plot_ss_choices = [60,floor(time_vector_p2d(end)*0.05),floor(time_vector_p2d(end)*0.5),time_vector_p2d(end)];
t_plot_ss_indices = find(ismember(time_vector_p2d,t_plot_ss_choices));

N = length(t_plot_ss_choices);
viridis_colors = viridis(N+1);
viridis_colors = viridis_colors(2:end,:);
for loop_idx = 1:length(t_plot_ss_choices)
    t_plot_idx = t_plot_ss_indices(loop_idx);
    ce_p2d_t_plot_idx = ce_p2d(t_plot_idx,:);
    ce_p2d_t_plot_idx_full_length = interp1(x_um,ce_p2d_t_plot_idx,x_um_full_length,'linear','extrap');
    plot(x_um_full_length,ce_p2d_t_plot_idx_full_length,'color',viridis_colors(loop_idx,:)); hold on;
end
annot_inset = floor(0.9*x_um_full_length(end));
text(annot_inset,1200,[num2str(t_plot_ss_choices(1)), ' s']);
text(annot_inset,1380,[num2str(t_plot_ss_choices(2)), ' s']);
text(annot_inset,1560,[num2str(t_plot_ss_choices(3)), ' s']);
text(annot_inset,1695,[num2str(t_plot_ss_choices(4)), ' s']);

%% 
line(1e6*[param_p2d{1}.len_p param_p2d{1}.len_p],get(gca,'YLim'),'linestyle','--','linewidth',1,'color','k');
line(1e6*[param_p2d{1}.len_p+param_p2d{1}.len_s param_p2d{1}.len_p+param_p2d{1}.len_s],get(gca,'YLim'),'linestyle','--','linewidth',1,'color','k');
% return;
%%
ylim([0 1750]);
no_of_y_tick_points = 8;
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
% box off;
MagInset(fig_h, -1, [1e6*(param_p2d{1}.len_p+param_p2d{1}.len_s*1.01), 1e6*(param_p2d{1}.len_p+param_p2d{1}.len_s+param_p2d{1}.len_n), 990 1080], [1e6*(param_p2d{1}.len_p+param_p2d{1}.len_s+param_p2d{1}.len_n*0.2), 1e6*(param_p2d{1}.len_p+param_p2d{1}.len_s+param_p2d{1}.len_n*0.8) 150 800], {'SW','NW';'SE','NE'});
ax_handle_inset = gca;
% ax_handle_inset.YAxis.TickValues = linspace(ax_handle_inset.YAxis.Limits(1),ax_handle_inset.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
set(ax_handle_inset,'xtick',[]);
set(ax_handle_inset,'ytick',[]);
d = findobj(gcf,'Type','text');
d(4).delete;
% allAxesInFigure = findall(fig_h,'type','axes');
annot_inset = floor(0.8*x_um_full_length(end));
text(annot_inset,1021,[num2str(t_plot_transient_choices(1)), ' s']);
text(annot_inset,1040,[num2str(t_plot_transient_choices(2)), ' s']);
text(annot_inset,1070,[num2str(t_plot_transient_choices(3)), ' s']);

shg;
% return;
%% Export to Tikz
extra_axis_options = 'yticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=2, /pgf/number format/fixed, }';
custom_m2t_fcn('ce_1C_at_various_times',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
close;
