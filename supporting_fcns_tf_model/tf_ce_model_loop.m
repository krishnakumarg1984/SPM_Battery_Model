%% Setup the parameters of the TF ce model
% load('Qen_tf_4p3z_scaled.mat');
% load('Qep_tf_4p3z_scaled.mat');
% load('Qes_tf_4p3z.mat');
% ce_tf_subsys_params = define_ce_tf_model_parameters(soc_init_pct);
extradata_for_solver.ce_tf_subsys_params= ce_tf_subsys_params;
extradata_for_solver.t0 = 0;
% Qen_tf_model = tf(Qen_tfest_scaled); clear Qen_tfest_scaled;
% Qep_tf_model = tf(Qep_tfest_scaled); clear Qep_tfest_scaled;
% Qes_tf_model = tf(Qes_tf_4p3z); clear Qes_tf_4p3z Qes_tfest_4p3z;
% t_end_common = min(t_finish,time_vector_p2d(end)-Ts);
extradata_for_solver.tf = t_finish;
% return;
%% Obtain constant parameters required for initialisation and simulation
% Retrieve the constant parameters of the model required for initialisation
% ce_init = ce_tf_subsys_params.ce_init;
% eps_n = ce_tf_subsys_params.eps_n;
% eps_s = ce_tf_subsys_params.eps_s;
% eps_p = ce_tf_subsys_params.eps_p;

% Ln = ce_tf_subsys_params.len_n;
% Ls = ce_tf_subsys_params.len_s;
% Lp = ce_tf_subsys_params.len_p;

%% Initialise the PP2D model's differential and algebraic variables at t = 0
% Qe_tf_vector_init = [eps_n*ce_init*Ln;eps_s*ce_init*Ls;eps_p*ce_init*Lp]; % at t = 0
a_tf_vector = [ce_init;... % a0
    0;...          % a2
    ce_init;...    % a3
    0;...          % a4
    0;...          % a5
    ce_init;...    % a6
    0]; % at t = 0 % a8

% num_iterations = ceil(t_end_common/Ts) + 1; % max no. of steps (assuming no cutoff)
% no_of_steps = ceil(t_finish/Ts) + 1; % max no. of steps (assuming no cutoff)
a_tf_vector_results = nan(length(a_tf_vector),no_of_steps);
a_tf_vector_results(:,1) = a_tf_vector;
% Qe_tf_vector_results = nan(length(Qe_tf_vector_init),no_of_steps);
% Qe_tf_vector_results(:,1) = Qe_tf_vector_init;

% tf_ce_sim_time_vector = nan(no_of_steps,1);
tf_ce_sim_time_vector = tf_interp_vector;
t0_local_tf = 0;
tf_local_tf = Ts;

% tf_ce_sim_time_vector(1) = t0_local_tf;
% load_current_vector      = nan(no_of_steps,1);
% sim_results_tf_ce = nan(no_of_steps,1);
% load_current_vector(1) = I_1C*(interp1(C_rate_profile(:,1),C_rate_profile(:,2),t0_local_tf,'previous','extrap'));
%% Run PP2D model simulation loop
progressbarText(0);
step_tf = 1;
while step_tf <= no_of_steps
%     I_load_t = I_1C*(interp1(C_rate_profile(:,1),C_rate_profile(:,2),tf_local_tf,'previous','extrap'));
    I_load_t = I_load_vector(step_tf);
%     load_current_vector(step_tf) = I_load_t;
    
    extradata_for_solver.load_curr_anonymousFcn = @(t_local,t0,tf) I_load_t; % This is an anonymous fcn (pure fcn)
    extradata_for_solver.a_vector = a_tf_vector;
    tspan_tf = [t0_local_tf tf_local_tf];
%     [~,Qes_vector_sol_tf] = ode45(@(t,y) rhs_tf_ce_sep_ODE(t,y,extradata_for_solver), tspan_tf, Qe_tf_vector_init(2));
%     Qen_vector_sol_tf = ((1/1000)*lsim(Qen_tf_model,I_load_t*ones(size(tspan_tf')),[0;1],[],'zoh')) + Qe_tf_vector_init(1);
%     Qes_vector_sol_tf = lsim(Qes_tf_model,I_load_t*ones(size(tspan_tf')),[0;1],[],'zoh') + Qe_tf_vector_init(2);
%     Qep_vector_sol_tf = ((1/1000)*lsim(Qep_tf_model,I_load_t*ones(size(tspan_tf')),[0;1],[],'zoh')) + Qe_tf_vector_init(3);
%     Qe_tf_vector_init = [Qen_vector_sol_tf(end);Qes_vector_sol_tf(end);Qep_vector_sol_tf(end)];

    [a_tf_vector,~,~] = solve_algebraic_a_vector(extradata_for_solver,Qe_tf_vector_results(:,step_tf));
    a_tf_vector_results(:,step_tf) = a_tf_vector;  % required for later post-processing
%     Qe_tf_vector_results(:,step_tf) = Qe_tf_vector_init; % Only for checking if an 8-th equation is possible or not
    
%     tf_ce_sim_time_vector(step_tf) = tf_local_tf;
    t0_local_tf = t0_local_tf + Ts;
    tf_local_tf = tf_local_tf + Ts;
    step_tf = step_tf + 1;
    progressbarText((step_tf-1)/no_of_steps);
end
extradata_for_solver = rmfield(extradata_for_solver,{'load_curr_anonymousFcn','a_vector'});
