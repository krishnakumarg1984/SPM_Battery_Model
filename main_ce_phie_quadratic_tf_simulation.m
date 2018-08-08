% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;
% warning('off','all');

%% Load user data, pre-process and run simulation loop
% load('p2d_sim_Jul_09_2018_14_04_28'); % 1C constant current dischg
% load('p2d_sim_Aug_07_2018_12_58_11'); % 1C constant current dischg with phie results
load('p2d_sim_Aug_07_2018_14_18_19'); % 1C constant current dischg with phie results with 30 nodes

% load('p2d_sim_Aug_07_2018_14_50_00'); % udds beginning at 50 percent soc with phie and 30 nodes each
run('user_inputs_for_sim.m');
run('pre_process_script.m');

run('quadratic_ce_model_loop.m');
run('postprocess_quadratic_results.m');

run('setup_tf_model_parameters.m');
run('tf_ce_model_loop.m');
run('postprocess_tf_results.m');
% return;
%%
run('compute_phie_op_tf.m');
% return;
%% Obtain ce interpolated at the right location
z_neg_newconvention = linspace(param_p2d{1}.len_n/(2*param_p2d{1}.Nn),param_p2d{1}.len_n - (param_p2d{1}.len_n/(2*param_p2d{1}.Nn)),param_p2d{1}.Nn);
z_sep_newconvention = linspace(param_p2d{1}.len_s/(2*param_p2d{1}.Ns),param_p2d{1}.len_s - (param_p2d{1}.len_s/(2*param_p2d{1}.Ns)),param_p2d{1}.Ns);
z_pos_newconvention = linspace(param_p2d{1}.len_p/(2*param_p2d{1}.Np),param_p2d{1}.len_p - (param_p2d{1}.len_p/(2*param_p2d{1}.Np)),param_p2d{1}.Np);

z_neg_p2d_plot = linspace(0,param_p2d{1}.len_n,Nn);
z_sep_p2d_plot = linspace(0,param_p2d{1}.len_s,Ns);
z_pos_p2d_plot = linspace(0,param_p2d{1}.len_p,Np);

x_cell_p2d_plot_um = [z_neg_p2d_plot,param_p2d{1}.len_n + z_sep_p2d_plot, param_p2d{1}.len_n + param_p2d{1}.len_s + z_pos_p2d_plot]*1e6; % um
x_cell_newconvention = [z_neg_newconvention, param_p2d{1}.len_n+z_sep_newconvention, param_p2d{1}.len_n + param_p2d{1}.len_s + z_pos_newconvention];
x_cell_newconvention_um = x_cell_newconvention*1e6;

ce_pos_p2d = ce_results_p2d(:,1:param_p2d{1}.Np);
ce_sep_p2d = ce_results_p2d(:,param_p2d{1}.Np+1:param_p2d{1}.Np+param_p2d{1}.Ns);
ce_neg_p2d = ce_results_p2d(:,param_p2d{1}.Np+param_p2d{1}.Ns+1:end);

ce_neg_p2d_newconvention = fliplr(ce_neg_p2d);
ce_sep_p2d_newconvention = fliplr(ce_sep_p2d);
ce_pos_p2d_newconvention = fliplr(ce_pos_p2d);

phie_pos_p2d = phie_results_p2d(:,1:param_p2d{1}.Np);
phie_sep_p2d = phie_results_p2d(:,param_p2d{1}.Np+1:param_p2d{1}.Np+param_p2d{1}.Ns);
phie_neg_p2d = phie_results_p2d(:,param_p2d{1}.Np+param_p2d{1}.Ns+1:end);

phie_diff_p2d_newconvention = phie_pos_p2d(:,1) - phie_neg_p2d(:,end);
% save_foldername = ['phie_op_results/', cellIdentifier, '/', load_profile_name];
% if exist(save_foldername,'dir')==0
%     mkdir(save_foldername);
% end
% save([save_foldername,'/phie_sysid_'...
%     datestr(now, 'mmm_dd_yyyy_HH_MM_SS')],'phie_op_vector_results','phie_op_term1_vector','phie_op_term2_vector'); % save workspace to file
% return;


%% Setup plots
golden_ratio = 1.618;
run('setup_line_colors.m'); % deletes axes/clears plots
cbrewerintergray = [189,189,189]/255;
cbrewerdarkgray = [99,99,99]/255;
cbrewer_Gnbubu_blue = [67,162,202]/255;
no_of_x_tick_points = 6;
no_of_y_tick_points = 6;

clc;
fig_h = clf;
fig_h.Units = 'centimeters';
fig_width_factor = 0.75; % scaling factor wrt to text width
figW_cm = 15.74776*fig_width_factor;     % textwidth (cm) reported by LaTeX doc with a scaling factor
figH_cm = figW_cm/golden_ratio;
if strcmp(load_profile_name,'cnst_dischg_soc_100_1C')
    set(0,'defaultaxesfontsize',12,'defaultlinelinewidth',2,'defaultpatchlinewidth',0.5,'defaultaxeslinewidth',1);
else
    set(0,'defaultaxesfontsize',12,'defaultlinelinewidth',1,'defaultpatchlinewidth',0.5,'defaultaxeslinewidth',1);
end
fig_h.Position = [fig_h.Position(1),fig_h.Position(2),figW_cm,figH_cm];
movegui('center');

% %% Plot
plot(time_vector_p2d,phie_diff_p2d_newconvention,'color',cbrewerdarkgray);
hold on;
plot(tf_ce_sim_time_vector,phie_op_vector_results,'color',cbrewer_Gnbubu_blue);
hold off;
ylabel('V $\quad$');
set(get(gca,'ylabel'),'rotation',0);
title('$\Delta\phi_\mathrm{e}(t) = \phi_\mathrm{e,poscc}(t) - \phi_\mathrm{e,negcc}(t)$');
% lgd = legend('P2d','SysID','location','northeast');
% legend boxoff;
% xlabel('time (s)');
text_x = 0.325;
text(text_x,0.725,'P2d model','Units','Normalized');
text(text_x,0.53,'SysID model','Units','Normalized');

ax_handle = gca;
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
InSet = get(ax_handle, 'TightInset');
if strcmp(load_profile_name,'udds_soc_50')
    ylim([-0.2 0.5]);
    set(ax_handle, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
    MagInset(fig_h,-1,[600 700 -0.09 0.15],[175 875 0.25 0.48],{'NW','SW';'NE','SE'});
end

% return;

%%
extra_axis_options = 'xticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=0,/pgf/number format/fixed,/pgf/number format/fixed zerofill,},yticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=2, /pgf/number format/fixed, /pgf/number format/fixed zerofill,}, ylabel absolute, ylabel style={rotate=-90}';
if strcmp(load_profile_name,'cnst_dischg_soc_100_1C')
    custom_m2t_fcn('phie_delta_cnst_1C',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
elseif strcmp(load_profile_name,'udds_soc_50')
    custom_m2t_fcn('phie_delta_udds',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
end
close;
