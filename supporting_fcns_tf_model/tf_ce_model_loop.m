%% Setup the parameters of the TF ce model
extradata_for_solver.ce_tf_subsys_params= ce_tf_subsys_params;
extradata_for_solver.t0 = 0;
extradata_for_solver.tf = t_finish;

%% Initialise the PP2D model's differential and algebraic variables at t = 0
a_tf_vector = [ce_init;... % a0
    0;...          % a2
    ce_init;...    % a3
    0;...          % a4
    0;...          % a5
    ce_init;...    % a6
    0]; % at t = 0 % a8

a_tf_vector_results = nan(length(a_tf_vector),no_of_steps);
a_tf_vector_results(:,1) = a_tf_vector;

tf_ce_sim_time_vector = tf_interp_vector;
t0_local_tf = 0;
tf_local_tf = Ts;

%% Run PP2D model simulation loop
progressbarText(0);
step_tf = 1;
while step_tf <= no_of_steps
    I_load_t = I_load_vector(step_tf);
    
    extradata_for_solver.load_curr_anonymousFcn = @(t_local,t0,tf) I_load_t; % This is an anonymous fcn (pure fcn)
    extradata_for_solver.a_vector = a_tf_vector;
    tspan_tf = [t0_local_tf tf_local_tf];

    [a_tf_vector,~,~] = solve_algebraic_a_vector(extradata_for_solver,Qe_tf_vector_results(:,step_tf));
    a_tf_vector_results(:,step_tf) = a_tf_vector;  % required for later post-processing
    
    t0_local_tf = t0_local_tf + Ts;
    tf_local_tf = tf_local_tf + Ts;
    step_tf = step_tf + 1;
    progressbarText((step_tf-1)/no_of_steps);
end
extradata_for_solver = rmfield(extradata_for_solver,{'load_curr_anonymousFcn','a_vector'});
