%% Set up code infrastructure required for LIONSIMBA (P2D) simulation (user has option to tune these parameters if desired. Otherwise use these defaults)
% P2D model
param_p2d{1} = Parameters_init_suppliedSOC_pct(soc_init_pct); % Refer to LIONSIMBA documentation for details of these settings
param_p2d{1}.SolidPhaseDiffusion    = 3;
param_p2d{1}.sim_datalog_interval   = Ts;
param_p2d{1}.Scope                  = 0;
param_p2d{1}.suppress_status_prints = 0;
param_p2d{1}.PrintHeaderInfo = 0;
param_p2d{1}.Nr_p = 15;param_p2d{1}.Nr_n = 15;

%% Initialise P2D model
initialState_p2d.Y  = [];
initialState_p2d.YP = [];
exit_reason_p2d     = 0; % This flag is used to notify the reason why a simulation came to a halt. If everything went well, it's value remains zero

%% Allocate storage & populate initial values of SPM state vector as well as result vectors
num_iterations = ceil(t_finish/Ts)+1;

%% Perform single-shot constant simulation
t_local_start  = 0;
progressbarText(0);
I_1C_p2d = param_p2d{1}.i_1C_density*param_p2d{1}.overall_surface_area_for_given_layers;

I_load = I_1C_p2d*C_rate;
fprintf('Simulating P2D model ... \n');

results_p2d = startSimulation(t_local_start,t_finish,[],-I_load/param_p2d{1}.overall_surface_area_for_given_layers,param_p2d);
load_curr_vector_p2d(:) = I_load;
time_vector_p2d = results_p2d.time{1};
cell_voltage_results_p2d = results_p2d.Voltage{1};
soc_pct_results_p2d = results_p2d.SOC{1};
if param_p2d{1}.SolidPhaseDiffusion==1 || param_p2d{1}.SolidPhaseDiffusion==2
    cs_avg_neg_results_p2d  = mean(results_p2d.cs_average{1}(:,param_p2d{1}.Np+1:end),2);
else
    cs_avg_neg_results_p2d  = mean(results_p2d.cs_average{1}(:,param_p2d{1}.Np*param_p2d{1}.Nr_p+1:end),2);
end
