x_cell = [z_pos,z_pos(end)+z_sep,z_pos(end)+z_sep(end)+z_neg];
last_index_to_plot = length(t_solved);

%% Setting up a few plot colors
color_imp_blue = [0,62,116]/255;
color_orange = [210,64,0]/255;
color_brick = [165,25,0]/255;
color_plum = [50,30,109]/255;
color_darkgreen = [2,137,59]/255;
color_tangerine = [236,115,0]/255;

%% Electrolyte Concentration over time at different spatial locations
cd('plot_codes_and_plots/');

close all;
figure; clf;
maximize(gcf);

load_curr_plotposition = [0.055 0.75 0.9 0.2];
h_current = subplot('Position',load_curr_plotposition);
% stairs(C_rate_profile(:,1),C_rate_profile(:,2)*I_1C,'o-','linewidth',2.5,'color',color_darkgreen,'MarkerEdgeColor',color_tangerine,'MarkerFaceColor',color_tangerine); % for short durations
plot(C_rate_profile(:,1),C_rate_profile(:,2)*I_1C,'linewidth',2.5,'color',color_darkgreen,'MarkerEdgeColor',color_tangerine,'MarkerFaceColor',color_tangerine); % for short durations
ylim([min(C_rate_profile(:,2)*I_1C)-5 max(C_rate_profile(:,2)*I_1C)+5]);
xlabel('$\mathrm{time (s)}$','interpreter','latex');ylabel('$I_\mathrm{Load}(t)$','interpreter','latex');
title('Load Current','interpreter','latex');
clear load_curr_plotposition;

ce_neg_sep_plotposition = [0.055 0.4 0.4 0.25];
h1 = subplot('Position',ce_neg_sep_plotposition);
plot(t_solved,ce_neg_sep_poly(1:last_index_to_plot));
hold on;
plot(t_solved,ce_neg_sep_ML(1:last_index_to_plot));
plot(time_p2d,ce_neg_sep_p2d);
ylim([min([ce_neg_sep_poly(:);ce_neg_sep_ML(:);ce_neg_sep_p2d(:)])-5 max([ce_neg_sep_poly(:);ce_neg_sep_ML(:);ce_neg_sep_p2d(:)])+5]);
title('Electrolyte Conc. at Neg -- Sep interface','interpreter','latex');
ylabel('$C_e (\mathrm{Ln},t)$','interpreter','latex');
legend('quadratic poly','symbolic data mining','P2D','location','best');
hold off;
clear ce_neg_sep_plotposition;

ce_neg_cc_plotposition = [0.55 0.4 0.4 0.25];
h2 = subplot('Position',ce_neg_cc_plotposition);
plot(t_solved,ce_neg_cc_poly(1:last_index_to_plot));
hold on;
plot(t_solved,ce_neg_cc_ML(1:last_index_to_plot));
plot(time_p2d,ce_neg_cc_p2d);
ylim([min([ce_neg_cc_poly(:);ce_neg_cc_ML(:);ce_neg_cc_p2d(:)])-5 max([ce_neg_cc_poly(:);ce_neg_cc_ML(:);ce_neg_cc_p2d(:)])+5]);
title('Electrolyte Conc. at Neg -- CC interface','interpreter','latex');
ylabel('$C_e (0,t)$','interpreter','latex');
legend('quadratic poly','symbolic data mining','P2D','location','best');
clear ce_neg_cc_plotposition;
hold off;
shg;

ce_pos_sep_plotposition = [0.055 0.075 0.4 0.25];
h3 = subplot('Position',ce_pos_sep_plotposition);
plot(t_solved,ce_pos_sep_poly(1:last_index_to_plot));
hold on;
plot(t_solved,ce_pos_sep_ML(1:last_index_to_plot));
plot(time_p2d,ce_pos_sep_p2d);
ylim([min([ce_pos_sep_poly(:);ce_pos_sep_ML(:);ce_pos_sep_p2d(:)])-5 max([ce_pos_sep_poly(:);ce_pos_sep_ML(:);ce_pos_sep_p2d(:)])+5]);
title('Electrolyte Conc. at Pos -- Sep interface','interpreter','latex');
xlabel('$\mathrm{time (s)}$','interpreter','latex');ylabel('$C_e (\mathrm{Ln + Ls},t)$','interpreter','latex');
legend('quadratic poly','symbolic data mining','P2D','location','best');
clear ce_pos_sep_plotposition;
hold off;

ce_pos_cc_plotposition = [0.55 0.075 0.4 0.25];
h4 = subplot('Position',ce_pos_cc_plotposition);
plot(t_solved,ce_pos_cc_poly(1:last_index_to_plot));
hold on;
plot(t_solved,ce_pos_cc_ML(1:last_index_to_plot));
plot(time_p2d,ce_pos_cc_p2d);
title('Electrolyte Conc. at Pos -- CC interface','interpreter','latex');
xlabel('$\mathrm{time (s)}$','interpreter','latex');ylabel('$C_e (\mathrm{Ln + Ls + Lp},t)$','interpreter','latex');
legend('quadratic poly','symbolic data mining','P2D','location','best');
clear ce_pos_cc_plotposition;
hold off;

linkaxes([h1,h2,h3,h4],'x');

set(gcf, 'Color', 'w');
% spatial_filename = ['temporal_ML_poly_p2d_ce_transient_' num2str(C_rate) 'C.png'];
specific_locs_over_time_filename = 'temporal_ML_poly_p2d_ce_transient_sim.png';
pause(2);
export_fig(specific_locs_over_time_filename);
clear specific_locs_over_time_filename
% close(gcf);

clc;
cd('../');