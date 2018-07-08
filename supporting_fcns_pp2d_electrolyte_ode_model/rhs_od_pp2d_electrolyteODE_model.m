function rhs_odes_pp2d_electrolyte = rhs_od_pp2d_electrolyteODE_model(t,y,extradata_for_solver)

% extradata_for_solver is used to extract the latest known a_vector components
% (Ordered) List of differential variables :
%  1. Qen
%  2. Qes
%  3. Qep

%% Retrieve the constant parameters of the model required by differential equations
tplus = extradata_for_solver.model_params.tplus;
A = extradata_for_solver.model_params.A;
F = extradata_for_solver.model_params.F;

Ln = extradata_for_solver.model_params.len_n;
Ls = extradata_for_solver.model_params.len_s;
Lp = extradata_for_solver.model_params.len_p;

De_eff_neg = extradata_for_solver.model_params.De_eff_neg;
De_eff_sep = extradata_for_solver.model_params.De_eff_sep;
De_eff_pos = extradata_for_solver.model_params.De_eff_pos;

%% Obtain the time-dependent 'parameters/coefficients' required to evaluate the RHS of each ode in the system
t0 = extradata_for_solver.t0;
tf = extradata_for_solver.tf;

I = extradata_for_solver.load_curr_anonymousFcn(t,t0,tf); % the RHS evaluates the applied load current at present solver-time 't'

a1 = extradata_for_solver.a_od_pp2d_vector(2);
a4 = extradata_for_solver.a_od_pp2d_vector(5);
a6 = extradata_for_solver.a_od_pp2d_vector(7);

%% RHS of individual ODE equations
rhs_od_pp2d_dQen_dt = ( (1 - tplus)*I/(A*F)) + 2*a1*De_eff_neg/Ln;
rhs_od_pp2d_dQes_dt = 2*a4*De_eff_sep/Ls;
rhs_od_pp2d_dQep_dt = (-(1 - tplus)*I/(A*F)) + 2*a6*De_eff_pos/Lp;

%% Assemble the RHS of the system
rhs_odes_pp2d_electrolyte = [rhs_od_pp2d_dQen_dt;...
                             rhs_od_pp2d_dQes_dt;...
                             rhs_od_pp2d_dQep_dt];
