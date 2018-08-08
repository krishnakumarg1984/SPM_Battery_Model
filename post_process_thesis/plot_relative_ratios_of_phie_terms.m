set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',1.5,'defaultlinelinewidth',2,'defaultpatchlinewidth',2);
golden_ratio = 1.618;

fig_width_factor = 0.75;
ax_width_factor = 0.725;
fig_ht_factor = 0.8;


load('phie_sysid_Aug_07_2018_22_38_15');
t_vector = linspace(0,length(phie_op_vector_results)-1,length(phie_op_vector_results));
run('setup_line_colors.m'); % deletes axes/clears plots
figW_cm = 15.74776*fig_width_factor;     % textwidth (cm) reported by LaTeX doc with a scaling factor
figH_cm = figW_cm/golden_ratio;
fig_h = clf;
fig_h.Units = 'centimeters';
fig_h.Position = [fig_h.Position(1),fig_h.Position(2),figW_cm,figH_cm];

cbrewerintergray = [189,189,189]/255;
cbrewerdarkgray = [99,99,99]/255;
cbrewer_Gnbubu_blue = [67,162,202]/255;
movegui('center');

%% Post-process and plot
no_of_x_tick_points = 9;
no_of_y_tick_points = 6;

clc;

plot(t_vector,phie_op_term1_vector,'color',cbrewerdarkgray); 
hold on;
plot(t_vector,phie_op_term2_vector,'color',cbrewerintergray);
plot(t_vector,phie_op_vector_results,'color','k');
hold off;
ylabel('V');
% lgd = legend('Concentration Polarization','Ohmic Loss in Bulk Solution','Overall $\Delta \phi_\mathrm{e}(t)$ in Electrolyte','location','northeast');
% legend boxoff;
text_x = 0.4675;
text(text_x,0.52,'Concentration Polarisation','Units','Normalized');
text(text_x,0.73,'Ohmic Loss in Bulk Solution','Units','Normalized');
text(text_x,0.21,'Overall $\Delta \phi_\mathrm{e}(t)$ in Electrolyte','Units','Normalized');

xlim([0 4000]);
ax_handle = gca;
ax_handle.YAxis.TickValues = linspace(ax_handle.YAxis.Limits(1),ax_handle.YAxis.Limits(2),no_of_y_tick_points); % not for voltage
ax_handle.XAxis.TickValues = linspace(ax_handle.XAxis.Limits(1),ax_handle.XAxis.Limits(2),no_of_x_tick_points);
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); % remove scientific multipliers in x-axis format
xlabel('time (s)');
ax_handle = gca;
InSet = get(ax_handle, 'TightInset');
set(ax_handle, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);

% return;

extra_axis_options = 'legend style={font=\footnotesize},title style={yshift=-1.75ex,},xticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=0,/pgf/number format/fixed,/pgf/number format/fixed zerofill,},yticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=2, /pgf/number format/fixed, }, ylabel absolute,  ylabel style={rotate=-90}';
custom_m2t_fcn('contribution_to_phie_1C',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
% close;

% return;
