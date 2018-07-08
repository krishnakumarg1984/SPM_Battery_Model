function [a_vector,A,b] = solve_algebraic_a_vector_od_poly(extradata_for_solver,Qe_vector_for_alg_system)

% extradata_for_solver is used only for extracting constant model parameters

% The components of the (ordered) algebraic vector are:
%  1. a0
%  2. a1
%  3. a2
%  4. a3
%  5. a4
%  6. a5
%  7. a6

%% Retrieve the constant parameters of the model required by algebraic equations
Ln = extradata_for_solver.model_params.len_n;
Ls = extradata_for_solver.model_params.len_s;
Lp = extradata_for_solver.model_params.len_p;

De_eff_neg = extradata_for_solver.model_params.De_eff_neg;
De_eff_sep = extradata_for_solver.model_params.De_eff_sep;
De_eff_pos = extradata_for_solver.model_params.De_eff_pos;

eps_n = extradata_for_solver.model_params.eps_n;
eps_s = extradata_for_solver.model_params.eps_s;
eps_p = extradata_for_solver.model_params.eps_p;

ce_init = extradata_for_solver.model_params.ce_init;

Qen_init = ce_init*eps_n*Ln;
Qes_init = ce_init*eps_s*Ls;
Qep_init = ce_init*eps_p*Lp;

%% Retrieve the components of the state vector from the latest solution
Qen = Qe_vector_for_alg_system(1);
Qes = Qe_vector_for_alg_system(2);
Qep = Qe_vector_for_alg_system(3);

%% formulate the system matrix
A = zeros(7,7);

%                a0               a1        a2              a3                a4        a5                a6
A(1,:) = [        1                1        -1               0                 0         0                 0 ];
A(2,:) = [        0                0         1               1                 1        -1                -1 ];
A(3,:) = [        0  2*De_eff_neg/Ln         0  -De_eff_sep/Ls                 0         0                 0 ];
A(4,:) = [        0                0         0   De_eff_sep/Ls   2*De_eff_sep/Ls         0   2*De_eff_pos/Lp ];
A(5,:) = [ eps_n*Ln     (eps_n*Ln)/3         0               0                 0         0                 0 ];
A(6,:) = [        0                0  eps_s*Ls    0.5*eps_s*Ls      (eps_s*Ls)/3         0                 0 ];
A(7,:) = [        0                0         0               0                 0  eps_p*Lp      (eps_p*Lp)/3 ];

b = [0;0;0;0;Qen-Qen_init;Qes-Qes_init;Qep-Qep_init];

a_vector = A\b;
end