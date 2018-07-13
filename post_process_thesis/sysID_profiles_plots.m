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

%% Random Binary Input Signal
N = 200;
u1 = idinput(N,'rbs',[],range);

%% Chirp
t_chirp_start = 0;
t_chirp_end = 3*N*Ts;
t = linspace(t_chirp_start,t_chirp_end,3*N);
f0 = 0;
f1 = 1e-1;
u2 = max(range)*chirp(t,f0,t(end),f1)';
% plot(t, u2, 'b*-');shg;
%% Periodic Random Binary Input Signal
bin_seq_Period = N; % seconds
bin_seq_Period_samples = ceil(bin_seq_Period/Ts); % seconds

bin_NumPeriod = 2;
u3 = idinput([bin_seq_Period_samples,NumChannel,bin_NumPeriod],'rbs',[],range);

u_train = [u1;u2;u3];
%% Training set is done. Plot setup
figW_cm = 15.74776*fig_width_factor;     % textwidth (cm) reported by LaTeX doc with a scaling factor
figH_cm = figW_cm/golden_ratio;
set(0,'defaultlinelinewidth',1.5,'defaultpatchlinewidth',0.5,'defaultaxeslinewidth',1);

N_excit_train = 3; % three types of excitation inputs
% divergingmap_plot = viridis(N+1);
divergingmap_plot = brewermap(N_excit_train+1,'PuOr');
% divergingmap_plot = divergingmap_plot(2:end,:);
% divergingmap_plot = [divergingmap_plot(1,:); divergingmap_plot(3:end,:)];
t_rbs = linspace(0,length(u1)-Ts,length(u1));
t_chirp = linspace(length(u1),length(u1)+length(u2)-Ts,length(u2));
t_perrbs = linspace(length(u1)+length(u2),length(u1)+length(u2)+length(u3)-Ts,length(u3));
% t_train = ([t_rbs t_chirp t_perrbs])';
%% Plot it now! 
fig_h = clf;
fig_h.Units = 'centimeters';
fig_h.Position = [fig_h.Position(1),fig_h.Position(2),figW_cm,figH_cm];
plot(t_rbs,u1,'color',divergingmap_plot(1,:),'linewidth',0.5); hold on;
plot(t_chirp,u2,'color',divergingmap_plot(3,:),'linewidth',1);
plot(t_perrbs,u3,'color',divergingmap_plot(2,:),'linewidth',0.5);hold off;
ylabel('$I_\mathrm{train}\, (\mathrm{A})$');
xlabel('time (sec)');
ax_handle = gca;
InSet = get(ax_handle, 'TightInset');
set(ax_handle, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)]);
annot_y = 0.95;
annot_fontsize = 14;
text(0.05,annot_y,'RBS','Units','normalized','FontSize',annot_fontsize);
text(0.4,annot_y,'Chirp','Units','normalized','FontSize',annot_fontsize);
text(0.725,annot_y,'Periodic RBS','Units','normalized','FontSize',annot_fontsize);
ylim([4*range(1),1.75*range(2)]);
yticks([-120,-80,-40,0,40,80,120]);
magoutside_range_factor = 1.1;
magplacement_ymin = -400;
magplacement_ymax = -200;
MagInset(fig_h, -1, [90, 140, magoutside_range_factor*range(1), magoutside_range_factor*range(2)], [120, 500,  magplacement_ymin, magplacement_ymax], {'SW','NW';'SE','NE'});
MagInset(fig_h, -1, [1000, 1050, magoutside_range_factor*range(1), magoutside_range_factor*range(2)], [700, 1080,  magplacement_ymin, magplacement_ymax], {'SW','NW';'SE','NE'});
annot_objs = findobj(gcf,'Type','text');
annot_objs(1).delete;
annot_objs(2).delete;
annot_objs(3).delete;
annot_objs(4).delete;
annot_objs(5).delete;
annot_objs(6).delete;
%% Export to Tikz
extra_axis_options = 'xticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=0,/pgf/number format/fixed,/pgf/number format/fixed zerofill,},yticklabel style={/pgf/number format/1000 sep=, /pgf/number format/precision=2, /pgf/number format/fixed, }';
return;
custom_m2t_fcn('sysid_train_input',[figW_cm,figH_cm]*10,[],false,extra_axis_options);
close;
return;
%% Random Guassian Signal
rgs_Period = N; % seconds
rgs_Period_samples = ceil(rgs_Period/Ts);
rgs_NumPeriod = 4;
u4 = idinput([rgs_Period_samples,NumChannel,rgs_NumPeriod],'rgs',[],0.5*range);

%% PRBS
prbs_Period = N; % seconds
prbs_Period_samples = ceil(prbs_Period/Ts);

u5 = idinput(N,'prbs',[],range);

%% Sum of sines
sine_samples_per_Period = 2*N;
sine_NumPeriod = 3;
[u6,freq] = idinput([sine_samples_per_Period 1 sine_NumPeriod],'sine',[],range);

%% Split into training and validation data
I_load_fit = [u1;u2;u3];
t_vector_fit = linspace(0,Ts*length(I_load_fit)-1,length(I_load_fit))';
csvwrite('sys_id_profile_fit.csv',[t_vector_fit I_load_fit],1,0);

I_load_validate = [u4;u5;u6];
t_vector_validate = linspace(0,Ts*length(I_load_validate)-1,length(I_load_validate))';
csvwrite('sys_id_profile_validate.csv',[t_vector_validate I_load_validate],1,0);

%% 
% plot(u3,'o-');
% ylim([-2 2]);
% shg;


I_load_fit = iddata([],I_load_fit,Ts);
plot(I_load_fit);

I_load_validate = iddata([],I_load_validate,Ts);
plot(I_load_validate);
shg;