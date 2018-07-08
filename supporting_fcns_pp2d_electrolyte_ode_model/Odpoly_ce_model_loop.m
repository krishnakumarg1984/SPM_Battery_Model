%% Setup the PP2D model and solver parameters
od_pp2d_model_params = define_pp2d_parameters(initial_soc_pct); % define the constant parameters of the model
params_to_remove = {'a_i','Rp_n','Rp_p','Dns','Dps','theta_max_pos','theta_max_neg','theta_min_pos','theta_min_neg','R','k_p','k_n','cs_maxp','cs_maxn','sig_eff','SolidPhaseDiffusion','TemperatureEnabled','T','Keff_pos','Keff_sep','Keff_neg','CutoffVoltage','CutoverVoltage','CutoffSOC','CutoverSOC'};
extradata_for_solver.pp2d_model_params = rmfield(od_pp2d_model_params,params_to_remove);

clear pp2d_model_params params_to_remove;

%% Obtain constant parameters required for initialisation and simulation
% Retrieve the constant parameters of the model required for initialisation
ce_init = extradata_for_solver.pp2d_model_params.ce_init;
eps_n = extradata_for_solver.pp2d_model_params.eps_n;
eps_s = extradata_for_solver.pp2d_model_params.eps_s;
eps_p = extradata_for_solver.pp2d_model_params.eps_p;

Ln = extradata_for_solver.pp2d_model_params.len_n;
Ls = extradata_for_solver.pp2d_model_params.len_s;
Lp = extradata_for_solver.pp2d_model_params.len_p;

I_1C = extradata_for_solver.pp2d_model_params.i_1C_density*extradata_for_solver.pp2d_model_params.A;

extradata_for_solver.t0 = t0;
extradata_for_solver.tf = tf;
%% Initialise the PP2D model's differential and algebraic variables at t = 0
Qe_od_poly_vector_init = [eps_n*ce_init*Ln;eps_s*ce_init*Ls;eps_p*ce_init*Lp]; % at t = 0
% a_od_poly_vector = [ce_init;0;ce_init;0;0;ce_init;0]; % at t = 0
a_od_poly_vector = zeros(7,1); % at t = 0

a_od_poly_vector_results = nan(length(a_od_poly_vector),no_of_steps);
a_od_poly_vector_results(:,1) = a_od_poly_vector;
Qe_od_poly_vector_results = nan(length(Qe_od_poly_vector_init),no_of_steps);
Qe_od_poly_vector_results(:,1) = Qe_od_poly_vector_init;

%% Run PP2D model simulation loop
progBar_od_poly = ProgressBar(no_of_steps,'Title', 'OD PP2D sim','UpdateRate',10); % intialise a progressbar
t0_local_od_poly = t0;
tf_local_od_poly = t0 + Ts;
step_od_poly = 1;
while step_od_poly <= no_of_steps
    I_load_t = I_1C*(interp1(C_rate_profile(:,1),C_rate_profile(:,2),t0_local_od_poly,'previous','extrap'));
    extradata_for_solver.load_curr_anonymousFcn = @(t0_local,t0,tf) I_load_t; % This is an anonymous fcn (pure fcn)
    extradata_for_solver.a_vector = a_od_poly_vector;
    tspan_od_poly = [t0_local_od_poly tf_local_od_poly];
    [~,Qe_vector_sol_od_poly] = ode45(@(t,y) rhs_pp2d_electrolyteODE_model(t,y,extradata_for_solver), tspan_od_poly, Qe_od_poly_vector_init);
    Qe_od_poly_vector_init = Qe_vector_sol_od_poly(end,:)';
    
    [a_od_poly_vector,~,~] = solve_algebraic_a_vector_od_poly(extradata_for_solver,Qe_od_poly_vector_init);
    a_od_poly_vector_results(:,step_od_poly+1) = a_od_poly_vector;  % required for later post-processing
    Qe_od_poly_vector_results(:,step_od_poly+1) = Qe_od_poly_vector_init; % Only for checking if an 8-th equation is possible or not
    
    t0_local_od_poly = t0_local_od_poly + Ts;
    tf_local_od_poly = tf_local_od_poly + Ts;
    step_od_poly = step_od_poly + 1;
    progBar_od_poly([], [], []);
end
progBar_od_poly.release();fprintf('\n');

extradata_for_solver = rmfield(extradata_for_solver,{'load_curr_anonymousFcn','a_vector'});