% x_cell_plot = [z_pos,z_pos(end)+z_sep,z_pos(end)+z_sep(end)+z_neg]*1e6;
p2d_node_locs_plot = p2d_node_locs*1e6;

last_index_to_plot = length(t_solved);

ylim_calcs = [p2d_results.ce{1}(:);ce_ML_tplot3(:);ce_ML_tplot4(:);ce_ML_tplot5(:);ce_ML_tplot6(:);ce_poly_tplot3(:);ce_poly_tplot4(:);ce_poly_tplot5(:);ce_poly_tplot6(:)];

%% Electrolyte Concentration (over full cell thickness) at different time snapshots
% cd('..');
cd('plot_codes_and_plots/');

close all;
figure; clf;
% maximize(gcf);
% pause(2);

load_curr_plotposition = [0.055 0.75 0.9 0.2];
h_current = subplot('Position',load_curr_plotposition);
% stairs(C_rate_profile(:,1),C_rate_profile(:,2)*I_1C,'o-','linewidth',2.5,'color',color_darkgreen,'MarkerEdgeColor',color_tangerine,'MarkerFaceColor',color_tangerine); % for short durations
plot(C_rate_profile(:,1),C_rate_profile(:,2)*I_1C,'linewidth',2.5,'color',color_darkgreen,'MarkerEdgeColor',color_tangerine,'MarkerFaceColor',color_tangerine); % for short durations
ylim([min(C_rate_profile(:,2)*I_1C)-5 max(C_rate_profile(:,2)*I_1C)+5]);
xlabel('$\mathrm{time (s)}$','interpreter','latex');ylabel('$I_\mathrm{Load}(t)$','interpreter','latex');
title('Load Current','interpreter','latex');
xlim([t0 tf]);
clear load_curr_plotposition;

ce_tplot3_plotposition = [0.055 0.4 0.4 0.25];
h1 = subplot('Position',ce_tplot3_plotposition);
plot(x_cell_plot,ce_poly_tplot3,'color',color_plum);
box on;hold on;
plot(x_cell_plot,ce_ML_tplot3,'color',color_brick);
plot(p2d_node_locs_plot,ce_tplot3_p2d,'color',color_tangerine);
legend('quadratic poly','symbolic data mining','P2D','location','best');
% xlim([min(min(x_cell_plot(:),p2d_node_locs_plot(:))) max(max(x_cell_plot,p2d_node_locs_plot))]);
ylabel('$C_e (x_{cell})$','interpreter','latex');
title(['Electrolyte Concentration over cell length at t = ' num2str(t_plot3) ' sec'],'interpreter','latex');
set(gca,'XTickLabel','');
clear ce_tplot3_plotposition;
% ylim([min(ylim_calcs)-5 max(ylim_calcs)+5]);

ce_tplot4_plotposition = [0.55 0.4 0.4 0.25];
h2 = subplot('Position',ce_tplot4_plotposition);
plot(x_cell_plot,ce_poly_tplot4,'color',color_plum);
box on;hold on;
plot(x_cell_plot,ce_ML_tplot4,'color',color_brick);
plot(p2d_node_locs_plot,ce_tplot4_p2d,'color',color_tangerine);
legend('quadratic poly','symbolic data mining','P2D','location','best');
% xlim([min(min(x_cell_plot,p2d_node_locs_plot)) max(max(x_cell_plot,p2d_node_locs_plot))]);
ylabel('$C_e (x_{cell})$','interpreter','latex');
title(['Electrolyte Concentration over cell length at t = ' num2str(t_plot4) ' sec'],'interpreter','latex');
set(gca,'XTickLabel','');
clear ce_tplot4_plotposition;
% ylim([min(ylim_calcs)-5 max(ylim_calcs)+5]);

ce_tplot5_plotposition = [0.055 0.075 0.4 0.25];
h3 = subplot('Position',ce_tplot5_plotposition);
plot(x_cell_plot,ce_poly_tplot5,'color',color_plum);box on;
hold on;
plot(x_cell_plot,ce_ML_tplot5,'color',color_brick);
plot(p2d_node_locs_plot,ce_tplot5_p2d,'color',color_tangerine);
legend('quadratic poly','symbolic data mining','P2D','location','best');
% xlim([min(min(x_cell_plot,p2d_node_locs_plot)) max(max(x_cell_plot,p2d_node_locs_plot))]);
xlabel('$x_{cell} (\mu \mathrm{m})$','interpreter','latex');
ylabel('$C_e (x_{cell})$','interpreter','latex');
title(['Electrolyte Concentration over cell length at t = ' num2str(t_plot5) ' sec'],'interpreter','latex');
clear ce_tplot5_plotposition;
% ylim([min(ylim_calcs)-5 max(ylim_calcs)+5]);

ce_tplot6_plotposition = [0.55 0.075 0.4 0.25];
h4 = subplot('Position',ce_tplot6_plotposition);
plot(x_cell_plot,ce_poly_tplot6,'color',color_plum);box on;
hold on;
plot(x_cell_plot,ce_ML_tplot6,'color',color_brick);
plot(p2d_node_locs_plot,ce_tplot6_p2d,'color',color_tangerine);
legend('quadratic poly','symbolic data mining','P2D','location','best');
% xlim([min(min(x_cell_plot,p2d_node_locs_plot)) max(max(x_cell_plot,p2d_node_locs_plot))]);
xlabel('$x_{cell} (\mu \mathrm{m})$','interpreter','latex');
ylabel('$C_e (x_{cell})$','interpreter','latex');
title(['Electrolyte Concentration over cell length at t = ' num2str(t_plot6) ' sec'],'interpreter','latex');
clear ce_tplot6_plotposition;
% ylim([min(ylim_calcs)-5 max(ylim_calcs)+5]);

linkaxes([h1 h2 h3 h4],'x');

%%
set(gcf, 'Color', 'w');
spatial_filename = ['spatial_ML_poly_p2d_ce_transient.png'];
pause(2);
export_fig(spatial_filename);
% close all;

clear p2d_node_locs_plot x_cell_plot;
clc;
cd('..');shg;return;