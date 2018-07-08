function [a_vector,A,b] = solve_algebraic_a_vector(extradata_for_solver,Qe_vector_for_alg_system)
% The components of the (ordered) algebraic vector are:
%  1. a0
%  2. a2
%  3. a3
%  4. a4
%  5. a5
%  6. a6
%  7. a8

%% Retrieve the constant parameters of the model required by algebraic equations
Ln = extradata_for_solver.ce_quadratic_subsys_params.len_n;
Ls = extradata_for_solver.ce_quadratic_subsys_params.len_s;
Lp = extradata_for_solver.ce_quadratic_subsys_params.len_p;

De_eff_neg = extradata_for_solver.ce_quadratic_subsys_params.De_eff_neg;
De_eff_sep = extradata_for_solver.ce_quadratic_subsys_params.De_eff_sep;
De_eff_pos = extradata_for_solver.ce_quadratic_subsys_params.De_eff_pos;

eps_n = extradata_for_solver.ce_quadratic_subsys_params.eps_n;
eps_s = extradata_for_solver.ce_quadratic_subsys_params.eps_s;
eps_p = extradata_for_solver.ce_quadratic_subsys_params.eps_p;

%% Retrieve the components of the state vector from the latest solution
Qen = Qe_vector_for_alg_system(1);
Qes = Qe_vector_for_alg_system(2);
Qep = Qe_vector_for_alg_system(3);

%% formulate the system matrix
A = nan(7,7);

%                a0               a2        a3              a4                a5        a6                a8
A(1,:) = [        1             Ln^2        -1               0                 0         0                 0 ];
A(2,:) = [        0                0        -1             -Ls             -Ls^2         1              Lp^2 ];
A(3,:) = [        0  2*Ln*De_eff_neg         0     -De_eff_sep                 0         0                 0 ];
A(4,:) = [        0                0         0     -De_eff_sep  -2*Ls*De_eff_sep         0  -2*Lp*De_eff_pos ];
A(5,:) = [ eps_n*Ln   (eps_n*Ln^3)/3         0               0                 0         0                 0 ];
A(6,:) = [        0                0  eps_s*Ls  0.5*eps_s*Ls^2    (eps_s*Ls^3)/3         0                 0 ];
A(7,:) = [        0                0         0               0                 0  eps_p*Lp    (eps_p*Lp^3)/3 ];


b = [0;0;0;0;Qen;Qes;Qep];

a_vector = A\b;
end