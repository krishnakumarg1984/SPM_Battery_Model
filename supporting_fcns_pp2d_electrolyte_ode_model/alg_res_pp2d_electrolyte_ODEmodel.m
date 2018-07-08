function res_alg_pp2d_electrolyte = alg_res_pp2d_electrolyte_ODEmodel(a_vector,extradata_for_solver,Qe_vector_for_alg_system)
% The components of the (ordered) algebraic state vector are:
%  1. a0
%  2. a1
%  3. a2
%  4. a3
%  5. a4
%  6. a5
%  7. a6

%% Retrieve the constant parameters of the model required by algebraic equations
Ln = extradata_for_solver.pp2d_model_params.len_n;
Ls = extradata_for_solver.pp2d_model_params.len_s;
Lp = extradata_for_solver.pp2d_model_params.len_p;

De_eff_neg = extradata_for_solver.pp2d_model_params.De_eff_neg;
De_eff_sep = extradata_for_solver.pp2d_model_params.De_eff_sep;
De_eff_pos = extradata_for_solver.pp2d_model_params.De_eff_pos;

eps_n = extradata_for_solver.pp2d_model_params.eps_n;
eps_s = extradata_for_solver.pp2d_model_params.eps_s;
eps_p = extradata_for_solver.pp2d_model_params.eps_p;

%% Retrieve the components of the state vector from the latest solution
Qen = Qe_vector_for_alg_system(1);
Qes = Qe_vector_for_alg_system(2);
Qep = Qe_vector_for_alg_system(3);

%% Separate the algebraic components from the initial guess 'a_vector'
a_vector = a_vector + randn(7,1);
a0 = a_vector(1);
a1 = a_vector(2);
a2 = a_vector(3);
a3 = a_vector(4);
a4 = a_vector(5);
a5 = a_vector(6);
a6 = a_vector(7);

%% Compute the residual of the indiviual algebraic equations within the DAE system
res_algEq1_electrolyte = a1*Ln^2 + a0 - a2;
res_algEq2_electrolyte = a6*Lp^2 + a5 - a4*Ls^2 - a3*Ls - a2;
res_algEq3_electrolyte = 2*a1*Ln*De_eff_neg - a3*De_eff_sep;
res_algEq4_electrolyte = -2*a6*Lp*De_eff_pos - (2*a4*Ls + a3)*De_eff_sep;
res_algEq5_electrolyte = Qen - (eps_n*((a1*Ln^3)/3 + a0*Ln));
res_algEq6_electrolyte = Qes - (eps_s*((a4*Ls^3)/3 + (a3*Ls^2/2) + a2*Ls));
res_algEq7_electrolyte = Qep - (eps_p*((a6*Lp^3)/3 + a5*Lp));

%% Assemble the individual residuals to form the residual array of the algebraic equations
res_alg_pp2d_electrolyte = [res_algEq1_electrolyte;...
                            res_algEq2_electrolyte;...
                            res_algEq3_electrolyte;...
                            res_algEq4_electrolyte;...
                            res_algEq5_electrolyte;...
                            res_algEq6_electrolyte;...
                            res_algEq7_electrolyte;...
                            ];
end