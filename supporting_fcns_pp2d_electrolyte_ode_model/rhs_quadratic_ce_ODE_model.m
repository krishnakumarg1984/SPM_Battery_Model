function rhs_odes_pp2d_electrolyte = rhs_quadratic_ce_ODE_model(t,y,extradata_for_solver)

% (Ordered) List of differential variables :
%  1. Qen
%  2. Qes
%  3. Qep

%% Retrieve the constant parameters of the model required by differential equations
tplus = extradata_for_solver.ce_quadratic_subsys_params.tplus;
A = extradata_for_solver.ce_quadratic_subsys_params.A;
F = extradata_for_solver.ce_quadratic_subsys_params.F;

Ln = extradata_for_solver.ce_quadratic_subsys_params.len_n;
Ls = extradata_for_solver.ce_quadratic_subsys_params.len_s;
Lp = extradata_for_solver.ce_quadratic_subsys_params.len_p;

De_eff_neg = extradata_for_solver.ce_quadratic_subsys_params.De_eff_neg;
De_eff_sep = extradata_for_solver.ce_quadratic_subsys_params.De_eff_sep;
De_eff_pos = extradata_for_solver.ce_quadratic_subsys_params.De_eff_pos;

%% Obtain the time-dependent 'parameters/coefficients' required to evaluate the RHS of each ode in the system
t0 = extradata_for_solver.t0;
tf = extradata_for_solver.tf;

I = extradata_for_solver.load_curr_anonymousFcn(t,t0,tf); % the RHS evaluates the applied load current at present solver-time 't'

a2 = extradata_for_solver.a_vector(2);
a5 = extradata_for_solver.a_vector(5);
a8 = extradata_for_solver.a_vector(7);

%% RHS of individual ODE equations
rhs_dQen_dt = ((1 - tplus)*I/(A*F)) + 2*a2*Ln*De_eff_neg;
rhs_dQes_dt = 2*a5*Ls*De_eff_sep;
rhs_dQep_dt = (-(1 - tplus)*I/(A*F)) + 2*a8*Lp*De_eff_pos;

%% Assemble the RHS of the system
rhs_odes_pp2d_electrolyte = [rhs_dQen_dt;rhs_dQes_dt;rhs_dQep_dt];
