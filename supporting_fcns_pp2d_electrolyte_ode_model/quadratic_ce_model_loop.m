%% Setup the parameters of the quadratic ce model
ce_quadratic_subsys_params = define_ce_quadratic_model_parameters(soc_init_pct);
extradata_for_solver.ce_quadratic_subsys_params= ce_quadratic_subsys_params;
extradata_for_solver.t0 = 0;
extradata_for_solver.tf = t_finish;
%% Obtain constant parameters required for initialisation and simulation
% Retrieve the constant parameters of the model required for initialisation
ce_init = ce_quadratic_subsys_params.ce_init;
eps_n = ce_quadratic_subsys_params.eps_n;
eps_s = ce_quadratic_subsys_params.eps_s;
eps_p = ce_quadratic_subsys_params.eps_p;

Ln = ce_quadratic_subsys_params.len_n;
Ls = ce_quadratic_subsys_params.len_s;
Lp = ce_quadratic_subsys_params.len_p;

%% Initialise the PP2D model's differential and algebraic variables at t = 0
Qe_quadratic_vector_init = [eps_n*ce_init*Ln;eps_s*ce_init*Ls;eps_p*ce_init*Lp]; % at t = 0
a_quadratic_vector = [ce_init;... % a0
    0;...          % a2
    ce_init;...    % a3
    0;...          % a4
    0;...          % a5
    ce_init;...    % a6
    0]; % at t = 0 % a8

num_iterations = ceil(t_finish/Ts) + 1; % max no. of steps (assuming no cutoff)
a_quadratic_vector_results = nan(length(a_quadratic_vector),num_iterations);
a_quadratic_vector_results(:,1) = a_quadratic_vector;
Qe_quadratic_vector_results = nan(length(Qe_quadratic_vector_init),num_iterations);
Qe_quadratic_vector_results(:,1) = Qe_quadratic_vector_init;

quadratic_ce_sim_time_vector = nan(num_iterations,1);
load_current_vector        = nan(num_iterations,1);
sim_results_quadratic_ce = nan(num_iterations,1);

%% Run PP2D model simulation loop
progressbarText(0);
t0_local_quadratic = 0;
tf_local_quadratic = Ts;
step_quadratic = 1;
while step_quadratic <= num_iterations
    I_load_t = I_1C*(interp1(C_rate_profile(:,1),C_rate_profile(:,2),t0_local_quadratic,'previous','extrap'));
    load_current_vector(step_quadratic) = I_load_t;
    
    extradata_for_solver.load_curr_anonymousFcn = @(t_local,t0,tf) I_load_t; % This is an anonymous fcn (pure fcn)
    extradata_for_solver.a_vector = a_quadratic_vector;
    tspan_quadratic = [t0_local_quadratic tf_local_quadratic];
    [~,Qe_vector_sol_quadratic] = ode45(@(t,y) rhs_quadratic_ce_ODE_model(t,y,extradata_for_solver), tspan_quadratic, Qe_quadratic_vector_init);
    Qe_quadratic_vector_init = Qe_vector_sol_quadratic(end,:)';
    
    [a_quadratic_vector,~,~] = solve_algebraic_a_vector(extradata_for_solver,Qe_quadratic_vector_init);
    a_quadratic_vector_results(:,step_quadratic+1) = a_quadratic_vector;  % required for later post-processing
    Qe_quadratic_vector_results(:,step_quadratic+1) = Qe_quadratic_vector_init; % Only for checking if an 8-th equation is possible or not
    
    t0_local_quadratic = t0_local_quadratic + Ts;
    tf_local_quadratic = tf_local_quadratic + Ts;
    step_quadratic = step_quadratic + 1;
    progressbarText((step_quadratic-1)/num_iterations);
end
extradata_for_solver = rmfield(extradata_for_solver,{'load_curr_anonymousFcn','a_vector'});
