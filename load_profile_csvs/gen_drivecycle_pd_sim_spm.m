clear; clc; close all;
set(0,'defaultaxesfontsize',12,'defaultaxeslinewidth',2,'defaultlinelinewidth',2.5,'defaultpatchlinewidth',2,'DefaultFigureWindowStyle','docked');
platform = 'BEV';
[Mv,Cd,Av,eta_drivetrain,no_of_cells_in_pack,overall_electrochemical_surface_area] = vehicleSpecs_for_drivecycle(platform);
V_init_equilibrium_truth = 3.8229; % volts (for Northrop cell)

%% UDDS
cycle_name = 'udds';

[P_density_vector_udds,P_cell_udds,~,~] = compute_drivecycle_power_density(cycle_name,Mv,Cd,Av,eta_drivetrain,no_of_cells_in_pack,overall_electrochemical_surface_area);

% P_density_vector = P_density_vector(1:30); % COMMENT THIS OUT LATER

%% 
cycle_name = 'eudc';
[P_density_vector_eudc,P_cell_eudc,~,~] = compute_drivecycle_power_density(cycle_name,Mv,Cd,Av,eta_drivetrain,no_of_cells_in_pack,overall_electrochemical_surface_area);

P_density_vector = [P_density_vector_udds;P_density_vector_eudc];

clear Cd Mv Av eta_drivetrain no_of_cells_in_pack overall_electrochemical_surface_area P_batt P_cell P_wheels;
max_NSteps_p2d  = length(P_density_vector);  % Max number of steps that may be required without considering voltage and SoC cutoffs and cutovers

param{1} = initialise_spm_model(V_init_equilibrium_truth);
Q = (5e-5*param{1}.cs_maxp)^2; % cannot simulate cts time system with process noise.
R = (2.5e-3)^2; % Variance of the measurement noise, v[k] (for voltage measurements)
sensor_noise = chol(R)'*randn(1); % sensor noise
sensor_noise = 0*sensor_noise;
V_init_equilibrium_meas = V_init_equilibrium_truth + sensor_noise;

%% Set up code infrastructure required for LIONSIMBA (P2D) simulation (user has option to tune these parameters if desired. Otherwise use these defaults)
Ts = 1; % All drivecycle data were captured at 1 sec intervals
param{1}.OperatingMode          = 2; % power density is the input
param{1}.SolidPhaseDiffusion    = 2;
param{1}.sim_datalog_interval   = Ts;
param{1}.Scope                  = 0;
param{1}.suppress_status_prints = 1; % if set to 0, might interfere with progressbar printing
param{1}.PrintHeaderInfo        = 0;
param{1}.Nr_p = 15;param{1}.Nr_n = 15;
initialState_p2d.Y  = [];
initialState_p2d.YP = [];

%% Perform an initial run to estimate overall number of runs needed (based on SoC drop at the end of run 1)
cycle_no = 1;
[param,~,~,~,~,local_SoCpctTrue_p2d_vector,~,~,~,~] = ...
    run_drive_cycle_loop(Ts,max_NSteps_p2d,cycle_no,P_density_vector,param,initialState_p2d,Q,R);
SoCpctTrue_init = 100*((param{1}.cs_p_init/param{1}.cs_maxp) - param{1}.theta_min_pos)/(param{1}.theta_max_pos - param{1}.theta_min_pos);  % True SOC percentage at initial time-step (simulation). For experimental data, this is known
if abs(local_SoCpctTrue_p2d_vector(end) - SoCpctTrue_init) < 1e-2
    fprintf('Too long simulation needed.... Aborting!\n');
    return;
else
    total_cycles_needed = ceil(SoCpctTrue_init/(abs(SoCpctTrue_init - local_SoCpctTrue_p2d_vector(end)))); % Conservative estimate of overall number of cycles_needed to deplete the cell
    % total_cycles_needed = 2; % COMMENT THIS OUT LATER
end
P_density_vector = P_density_vector(:, ones(total_cycles_needed,1));
P_density_vector = P_density_vector(:);
max_NSteps_p2d = total_cycles_needed*max_NSteps_p2d;
%%
% t_local_start_vector = NaN(max_NSteps_p2d*total_cycles_needed,1);
% t_local_start_vector(1) = 0;
% t_local_start_vector(2:local_idx) = local_t_local_start_vector(2:local_idx);
% cell_current_p2d_vector = NaN(max_NSteps_p2d*total_cycles_needed,1);
% cell_current_p2d_vector(1:local_idx) = local_cell_current_p2d_vector(1:local_idx);
% cs_avg_pos_p2d_vector = NaN(max_NSteps_p2d*total_cycles_needed,1);
% cs_avg_pos_p2d_vector(1) = param{1}.cs_p_init;% + chol(Q)'*randn(1);  % process noise added only in simulation, (since we don't know the true model?).  from t = -infinity to before simulation
% cs_avg_pos_p2d_vector(2:local_idx) = local_cs_avg_pos_p2d_vector(2:local_idx);
% SoCpctTrue_p2d_vector = NaN(max_NSteps_p2d*total_cycles_needed,1); % True SOC percentage for all time-steps (simulation). For experimental data, this is known
% SoCpctTrue_p2d_vector(1) = SoCpctTrue_init;
% SoCpctTrue_p2d_vector(2:local_idx) = local_SoCpctTrue_p2d_vector(2:local_idx);
% v_cell_truth_p2d_vector = NaN(max_NSteps_p2d*total_cycles_needed,1);
% v_cell_truth_p2d_vector(1:local_idx) = local_v_cell_truth_p2d_vector(1:local_idx);
% v_cell_meas_p2d_vector = NaN(max_NSteps_p2d*total_cycles_needed,1);
% v_cell_meas_p2d_vector(1:local_idx) = local_v_cell_meas_p2d_vector(1:local_idx);

[param,idx,t_local_start_vector,cell_current_p2d_vector,cs_avg_pos_p2d_vector,SoCpctTrue_p2d_vector,v_cell_truth_p2d_vector,v_cell_meas_p2d_vector,initialState_p2d,~] = ...
    run_drive_cycle_loop(Ts,max_NSteps_p2d,cycle_no,P_density_vector,param,initialState_p2d,Q,R);

%% For "true measurement data", Perform P2D simulation in a loop & log results every 'Ts' seconds
param{1}.JacobianFunction = [];
results_p2d.JacobianFun = [];
results_p2d.original = [];
results_p2d.parameters{1,1}.JacobianFunction = [];
t_local_start_vector = t_local_start_vector(1:idx);
cell_current_p2d_vector = cell_current_p2d_vector(1:idx);
cs_avg_pos_p2d_vector = cs_avg_pos_p2d_vector(1:idx);
SoCpctTrue_p2d_vector = SoCpctTrue_p2d_vector(1:idx);
v_cell_truth_p2d_vector = v_cell_truth_p2d_vector(1:idx);
v_cell_meas_p2d_vector = v_cell_meas_p2d_vector(1:idx);
time_p2d_vector = t_local_start_vector;
clear process_noise sensor_noise max_NSteps idx cs_avg_pos_p2d_intermediate t_local_start_vector;
clc;

%% Save Results to disk
clear V_init_equilibrium_truth;
clear cycle_name cycle_no local_SoCpctTrue_p2d_vector initialState_p2d P_cell_eudc P_cell_udds P_density_vector_eudc P_density_vector_udds;
save_filename = 'p2d_simulated_truth_drivecycle_udds_eudc.mat';
save(save_filename);
% dlmwrite('p2d_simulated_truth_drivecycle_udds_eudc.txt',[time_p2d_vector SoCpctTrue_p2d_vector P_density_vector cell_current_p2d_vector v_cell_meas_p2d_vector],'delimiter','\t');

%% Post-process & plot
close all;
h0 = subplot(411);
plot(time_p2d_vector,P_density_vector,'-');
ylabel('Appl p (W/m^2)');
h1 = subplot(412);
plot(time_p2d_vector,cell_current_p2d_vector,'-');
ylabel('Current [A]');
ylim([min(cell_current_p2d_vector)-5 max(cell_current_p2d_vector)+5]);
h2 = subplot(413);
plot(time_p2d_vector,SoCpctTrue_p2d_vector);
% hold on;
% plot(time,mean(SoCpctTrue)*ones(length(time),1));
% ylim([mean(SoCpctTrue)-1 mean(SoCpctTrue)+1]);
ylabel('SOC [%]');
h3 = subplot(414);
plot(time_p2d_vector,v_cell_meas_p2d_vector);
ylabel('Measured V');

linkaxes([h1 h2 h3],'x');
xlim([time_p2d_vector(1) time_p2d_vector(end)]);
linkaxes([h0 h1 h2 h3],'x');
clear h0 h1 h2 h3;
xlabel('Time [sec]');

figure(1);shg;



function [param,local_idx,local_t_local_start_vector,local_cell_current_p2d_vector,local_cs_avg_pos_p2d_vector,local_SoCpctTrue_p2d_vector,local_v_cell_truth_p2d_vector,local_v_cell_meas_p2d_vector,initialState_p2d,cycle_no] = run_drive_cycle_loop(Ts,max_NSteps_p2d,cycle_no,P_density_vector,param,initialState_p2d,Q,R)

%% Allocate storage for truth results
local_t_local_start_vector = NaN(max_NSteps_p2d,1);
local_cell_current_p2d_vector = NaN(max_NSteps_p2d,1);
local_cs_avg_pos_p2d_vector = NaN(max_NSteps_p2d+1,1); % One state variable only for now (cs_avg_pos)

local_SoCpctTrue_p2d_vector = NaN(max_NSteps_p2d+1,1); % True SOC percentage for all time-steps (simulation). For experimental data, this is known

local_v_cell_truth_p2d_vector = NaN(max_NSteps_p2d+1,1);
local_v_cell_meas_p2d_vector = NaN(max_NSteps_p2d+1,1);

progBar = ProgressBar(max_NSteps_p2d,'Title', 'Simulation running','UpdateRate',10); % intialise a progressbar

for local_idx = 1:max_NSteps_p2d  % Forward looking interval,i.e. start at k and propel forward. Results are logged at t_finish except at first time-step
    local_t_local_start_vector(local_idx) = (local_idx-1)*Ts;
    t_local_finish_point = local_t_local_start_vector(local_idx) + Ts;
    results_p2d = startSimulation(local_t_local_start_vector(local_idx),t_local_finish_point,initialState_p2d,-P_density_vector(local_idx),param);
    local_cell_current_p2d_vector(local_idx) = results_p2d.curr_density(1);
    local_SoCpctTrue_p2d_vector(local_idx+1) = results_p2d.SOC{1}(end);
    
    if param{1}.SolidPhaseDiffusion==1 || param{1}.SolidPhaseDiffusion==2
        cs_avg_pos_p2d_intermediate = mean(results_p2d.cs_average{1}(:,1:param{1}.Np),2);
        local_cs_avg_pos_p2d_vector(local_idx+1) = cs_avg_pos_p2d_intermediate(end);
    else
        cs_avg_pos_p2d_intermediate = mean(results_p2d.cs_average{1}(:,1:param{1}.Np*param{1}.Nr_p),2);
        local_cs_avg_pos_p2d_vector(local_idx+1) = cs_avg_pos_p2d_intermediate(end);
    end
    
    if local_idx==1 && cycle_no==1
        param{1}.JacobianFunction  = results_p2d.JacobianFun;
        local_v_cell_truth_p2d_vector(local_idx:local_idx+1) = [results_p2d.Voltage{1}(1) results_p2d.Voltage{1}(end)];
        sensor_noise = chol(R)'*randn(2,1); % sensor noise
        sensor_noise = 0*sensor_noise;
        local_v_cell_meas_p2d_vector(local_idx:local_idx+1) = local_v_cell_truth_p2d_vector(local_idx:local_idx+1) + sensor_noise;
    else
        local_v_cell_truth_p2d_vector(local_idx+1) = results_p2d.Voltage{1}(end); % output also seems to be getting updated at the next time-step, instead of this present time-step
        sensor_noise = chol(R)'*randn(1); % sensor noise
        sensor_noise = 0*sensor_noise;
        local_v_cell_meas_p2d_vector(local_idx+1) = local_v_cell_truth_p2d_vector(local_idx+1) + sensor_noise;
    end
    
    initialState_p2d = results_p2d.initialState;
    
    progBar([], [], []);
end

progBar.release();fprintf('\n');
clear progBar t_local_finish_point;
% cycle_no = cycle_no + 1;
end
% vim: set nospell nowrap textwidth=0 wrapmargin=0 formatoptions-=t