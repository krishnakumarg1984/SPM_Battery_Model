% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;

C_rate = 1;
run('user_inputs_for_sim_cnst_curr_chg');
run('single_shot_p2d_cnst.m'); % p2d simulation loop is in this script

%% Set-up x-grid for plot
x_pos = linspace(param_p2d{1}.len_p/(2*param_p2d{1}.Np), param_p2d{1}.len_p-(param_p2d{1}.len_p/(2*param_p2d{1}.Np)), param_p2d{1}.Np);
x_sep = param_p2d{1}.len_p + linspace(param_p2d{1}.len_s/(2*param_p2d{1}.Ns), param_p2d{1}.len_s-(param_p2d{1}.len_s/(2*param_p2d{1}.Ns)), param_p2d{1}.Ns);
x_neg = param_p2d{1}.len_p + param_p2d{1}.len_s + linspace(param_p2d{1}.len_n/(2*param_p2d{1}.Nn), param_p2d{1}.len_n-(param_p2d{1}.len_n/(2*param_p2d{1}.Nn)), param_p2d{1}.Nn);
x = [x_pos, x_sep, x_neg];
x_um = x*1e6;
x_um_full_length = 1e6*linspace(0,param_p2d{1}.len_p+param_p2d{1}.len_s+param_p2d{1}.len_n,2*(param_p2d{1}.Np+param_p2d{1}.Ns+param_p2d{1}.Nn));
%% transient plots
t_plot_transient_choices = [1,10,20];
t_plot_transient_indices = find(ismember(time_vector_p2d,t_plot_transient_choices));

clf;
N = length(t_plot_transient_choices);
OrRdmap = brewermap(N+1,'OrRd'); % OrRd upto 4 colors is a) colorblind-safe, b) print-safe and c)photocopy-safe
set(0,'DefaultAxesColorOrder',OrRdmap(2:end,:)); 
% set(0,'DefaultAxesColorOrder',viridis(N)); 
for loop_idx = 1:length(t_plot_transient_choices)
    t_plot_idx = t_plot_transient_indices(loop_idx);
    ce_p2d_t_plot_idx = ce_p2d(t_plot_idx,:);
    ce_p2d_t_plot_idx_full_length = interp1(x_um,ce_p2d_t_plot_idx,x_um_full_length,'linear','extrap');
    plot(x_um_full_length,ce_p2d_t_plot_idx_full_length); hold on;
end
xlim([0 x_um_full_length(end)]);

annot_x = floor(0.5*x_um_full_length(end));
text(annot_x,990,[num2str(t_plot_transient_choices(1)), ' s']);
text(annot_x,1290,[num2str(t_plot_transient_choices(2)), ' s']);
text(annot_x,1470,[num2str(t_plot_transient_choices(3)), ' s']);
% text(annot_x,1670,[num2str(t_plot_transient_choices(4)), ' s']);

shg;

% %% SS Plot
% t_plot_ss_choices = [60,floor(time_vector_p2d(end)*0.05),floor(time_vector_p2d(end)*0.5),time_vector_p2d(end)];
% t_plot_ss_indices = find(ismember(time_vector_p2d,t_plot_ss_choices));
% 
% % run('setup_line_colors.m'); % deletes axes/clears plots
% 
% clf;
% N = length(t_plot_ss_choices);
% % set(0,'DefaultAxesColorOrder',brewermap(N,'OrRd')); 
% set(0,'DefaultAxesColorOrder',viridis(N)); 
% for loop_idx = 1:length(t_plot_ss_choices)
%     t_plot_idx = t_plot_ss_indices(loop_idx);
%     ce_p2d_t_plot_idx = ce_p2d(t_plot_idx,:);
%     ce_p2d_t_plot_idx_full_length = interp1(x_um,ce_p2d_t_plot_idx,x_um_full_length,'linear','extrap');
%     plot(x_um_full_length,ce_p2d_t_plot_idx_full_length); hold on;
% end
% xlim([0 x_um_full_length(end)]);
% 
% annot_x = floor(0.9*x_um_full_length(end));
% text(annot_x,1100,[num2str(t_plot_ss_choices(1)), ' s']);
% text(annot_x,1290,[num2str(t_plot_ss_choices(2)), ' s']);
% text(annot_x,1470,[num2str(t_plot_ss_choices(3)), ' s']);
% text(annot_x,1670,[num2str(t_plot_ss_choices(4)), ' s']);
% 
% shg;