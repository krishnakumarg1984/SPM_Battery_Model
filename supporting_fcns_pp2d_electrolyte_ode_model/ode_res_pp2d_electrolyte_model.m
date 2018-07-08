function [res_odes_pp2d_electrolyte,rhs_odes_pp2d_electrolyte] = ode_res_pp2d_electrolyte_model(t,XZ,dXZ_dt,ida_user_data)

% (Ordered) List of differential variables :
%  1. Qen
%  2. Qes
%  3. Qep

Z = XZ(ida_user_data.n_diff+1:end); % Retrieve the purely-algebraic variables

%% Retrieve the constant parameters of the model required by differential equations
tplus = ida_user_data.pp2d_model_params.tplus;
A = ida_user_data.pp2d_model_params.A;
F = ida_user_data.pp2d_model_params.F;

Ln = ida_user_data.pp2d_model_params.len_n;
Ls = ida_user_data.pp2d_model_params.len_s;
Lp = ida_user_data.pp2d_model_params.len_p;

De_eff_neg = ida_user_data.pp2d_model_params.De_eff_neg;
De_eff_sep = ida_user_data.pp2d_model_params.De_eff_sep;
De_eff_pos = ida_user_data.pp2d_model_params.De_eff_pos;

%% Retrieve the required time-derivatives of state variables
dQen_dt = dXZ_dt(1);
dQes_dt = dXZ_dt(2);
dQep_dt = dXZ_dt(3);

%% Retrieve the required algebraic variables (components from the algebraic vector)
a1 = Z(2);
a4 = Z(5);
a6 = Z(7);

I = Z(8);

%% RHS of ODE equations
rhs_dQen_dt = ((1 - tplus)*I/(A*F)) + 2*a1*Ln*De_eff_neg;
rhs_dQes_dt = 2*a4*Ls*De_eff_sep;
rhs_dQep_dt = (-(1 - tplus)*I/(A*F)) + 2*a6*Lp*De_eff_pos;

rhs_odes_pp2d_electrolyte = [rhs_dQen_dt;rhs_dQes_dt;rhs_dQep_dt];

%% Compute the residual of the indiviual ODEs within the DAE system
res_dQen_dt = dQen_dt - rhs_dQen_dt;
res_dQes_dt = dQes_dt - rhs_dQes_dt;
res_dQep_dt = dQep_dt - rhs_dQep_dt;

%% Assemble the individual residuals to form the residual array of the ODE equations
res_odes_pp2d_electrolyte = [res_dQen_dt;res_dQes_dt;res_dQep_dt];
