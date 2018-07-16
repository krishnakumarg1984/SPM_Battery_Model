% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

clear;clc; format short g; format compact; close all;
fig_width_factor = 1; % scaling factor wrt to text width
golden_ratio = 1.618;

warning('off','Ident:dataprocess:idinput7');
soc_init_pct = 50;
param{1} = Parameters_init_suppliedSOC_pct(soc_init_pct);

I_1C = param{1}.overall_surface_area_for_given_layers*param{1}.I1C; % Amps
range = 2*[-I_1C I_1C];
NumChannel = 1;
Ts = 1;
N = 200;

%% Validation data set: Random Guassian Signal
rgs_Period = N; % seconds
rgs_Period_samples = ceil(rgs_Period/Ts);
rgs_NumPeriod = 4;
u4 = idinput([rgs_Period_samples,NumChannel,rgs_NumPeriod],'rgs',[],0.5*range);
u4(u4>range(2)) = range(2);
u4(u4<range(1)) = range(1);

%% PRBS
prbs_Period = N; % seconds
prbs_Period_samples = ceil(prbs_Period/Ts);

u5 = idinput(N,'prbs',[],range);

%% Sum of sines
sine_samples_per_Period = 2*N;
sine_NumPeriod = 3;
[u6,freq] = idinput([sine_samples_per_Period 1 sine_NumPeriod],'sine',[],range);

%% Validation dataset is done. Plot setup
figW_cm = 15.74776*fig_width_factor;     % textwidth (cm) reported by LaTeX doc with a scaling factor
figH_cm = figW_cm/golden_ratio;
set(0,'defaultlinelinewidth',1.5,'defaultpatchlinewidth',0.5,'defaultaxeslinewidth',1);

N_excit_validation = 3; % three types of validation inputs
% divergingmap_plot = viridis(N+1);
divergingmap_plot = brewermap(N_excit_validation+1,'PuOr');
% divergingmap_plot = divergingmap_plot(2:end,:);
% divergingmap_plot = [divergingmap_plot(1,:); divergingmap_plot(3:end,:)];
t_rgs = linspace(0,length(u4)-Ts,length(u4));
t_chirp = linspace(length(u4),length(u4)+length(u5)-Ts,length(u5));
t_perrbs = linspace(length(u4)+length(u5),length(u4)+length(u5)+length(u6)-Ts,length(u6));
% t_validation = ([t_rgs t_chirp t_perrbs])';
%% Plot it now! 
fig_h = clf;
fig_h.Units = 'centimeters';
fig_h.Position = [fig_h.Position(1),fig_h.Position(2),figW_cm,figH_cm];
plot(t_rgs,u4,'color',divergingmap_plot(1,:),'linewidth',0.5); hold on;
plot(t_chirp,u5,'color',divergingmap_plot(4,:),'linewidth',1);
plot(t_perrbs,u6,'color',divergingmap_plot(2,:),'linewidth',0.5);hold off;
ylabel('$I_\mathrm{train}\, (\mathrm{A})$');
xlabel('time (sec)');
ax_handle = gca;
InSet = get(ax_handle, 'TightInset');
set(ax_handle, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
annot_y = 0.925;
annot_fontsize = 14;
text(0.15,annot_y,'RGS','Units','normalized','FontSize',annot_fontsize);
text(0.3625,annot_y,'PRBS','Units','normalized','FontSize',annot_fontsize);
text(0.65,annot_y,'Multisine','Units','normalized','FontSize',annot_fontsize);
xlim([0 2200]);
ylim([4*range(1),1.75*range(2)]);
yticks([-120,-80,-40,0,40,80,120]);
magoutside_range_factor = 1.1;
magplacement_ymin = -400;
magplacement_ymax = -200;
MagInset(fig_h, -1, [350, 450, magoutside_range_factor*range(1), magoutside_range_factor*range(2)], [170, 640,  magplacement_ymin, magplacement_ymax], {'SW','NW';'SE','NE'});
MagInset(fig_h, -1, [930, 980, magoutside_range_factor*range(1), magoutside_range_factor*range(2)], [870, 1340,  magplacement_ymin, magplacement_ymax], {'SW','NW';'SE','NE'});
MagInset(fig_h, -1, [1675, 1775, magoutside_range_factor*range(1), magoutside_range_factor*range(2)], [1600, 2070,  magplacement_ymin, magplacement_ymax], {'SW','NW';'SE','NE'});
annot_objs = findobj(gcf,'Type','text');
annot_objs(1).delete;
annot_objs(2).delete;
annot_objs(3).delete;
annot_objs(4).delete;
annot_objs(5).delete;
annot_objs(6).delete;
annot_objs(7).delete;
annot_objs(8).delete;
annot_objs(9).delete;
%% Export to Tikz
extra_axis_options = 'xticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=0,/pgf/number format/fixed,/pgf/number format/fixed zerofill,},yticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=2, /pgf/number format/fixed, },ylabel style={xshift=1.75cm}';
return;
custom_m2t_fcn('sysid_validation_input',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
close;
return;

%% Split into training and validation data
I_load_fit = [u4;u5;u6];
t_vector_fit = linspace(0,Ts*length(I_load_fit)-1,length(I_load_fit))';
csvwrite('sys_id_profile_fit.csv',[t_vector_fit I_load_fit],1,0);

I_load_validate = [u4;u5;u6];
t_vector_validate = linspace(0,Ts*length(I_load_validate)-1,length(I_load_validate))';
csvwrite('sys_id_profile_validate.csv',[t_vector_validate I_load_validate],1,0);

%% 
% plot(u6,'o-');
% ylim([-2 2]);
% shg;


I_load_fit = iddata([],I_load_fit,Ts);
plot(I_load_fit);

I_load_validate = iddata([],I_load_validate,Ts);
plot(I_load_validate);
shg;