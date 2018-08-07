% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;
golden_ratio = 1.618;
ax_width_factor = 0.75;
fig_ht_factor = 0.75;

fig_width_factor = 1;
figW_cm = 15.74776*fig_width_factor; % textwidth (cm) reported by thesis latex doc from layouts package. Additonally include a scaling factor
figH_cm = figW_cm/golden_ratio;

clf;
load('p2d_sim_Aug_07_2018_14_18_19'); % cnst curr p2d sim with 30 nodes;
plot(time_vector_p2d,cell_voltage_results_p2d);
hold on;
load('cts_sim_May_26_2018_20_35_34'); % cnst curr 1C basic spm
plot(spm_sim_time_vector,v_cell_sim_results_spm);
load('cts_sim_Aug_07_2018_16_20_26'); % cnst curr 1C spm with phie
plot(spm_sim_time_vector,v_cell_sim_results_spm);
hold off;


% axes(ha(2*idx-1));
% xlabel('time (s)');
% axes(ha(2*idx));
% xlabel('time (s)');
return;
%% Export to Tikz
extra_axis_options = 'xticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=0,/pgf/number format/fixed,/pgf/number format/fixed zerofill,},yticklabel style={/pgf/number format/1000 sep=,},';
custom_m2t_fcn('const_curr_dischg_soc',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
% custom_m2t_fcn('const_curr_dischg_voltage',[figW_cm,figH_cm]*10,[],false,extra_axis_options);

close;
% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t: