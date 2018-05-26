function v_cell = spm_three_states_battery_voltage(x,u,param)
% returns v_cell given the vector x and input u at any given time-step
% x(1) = q_pos, x(2) = q_neg, x(3) = cs_avg_neg
% u = load current (A). Positive implies discharge.

% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

R      = param.R;
T      = param.Tref;
F      = param.F;

R_pos  = param.R_p;
R_neg  = param.R_n;

Ds_pos = param.D_p;
Ds_neg = param.D_n;

a_pos  = param.a_p;
a_neg  = param.a_n;

L_pos  = param.len_p;
L_neg  = param.len_n;

k_p    = param.k_p;
k_n    = param.k_n;

ce     = param.ce_init;
A      = param.A;

cs_max_n = param.cs_max_n;
cs_max_p = param.cs_max_p;

theta_min_neg = param.theta_min_neg;
theta_min_pos = param.theta_min_pos;
theta_max_neg = param.theta_max_neg;
theta_max_pos = param.theta_max_pos;

% Extract permissible bounds for solid concentrations
cs_surf_pos_lb = param.theta_max_pos*param.cs_max_p;
cs_surf_pos_ub = param.theta_min_pos*param.cs_max_p;
cs_surf_neg_lb = param.theta_min_neg*param.cs_max_n;
cs_surf_neg_ub = param.theta_max_neg*param.cs_max_n;

% Extract the function handles for Uocp_pos and Uocp_neg
compute_Uocp_pos = param.compute_Uocp_pos;
compute_Uocp_neg = param.compute_Uocp_neg;

%% Compute surface concentrations
cs_surf_neg = x(3) + (8*R_neg/35)*x(2) - (R_neg/(35*Ds_neg*a_neg*L_neg*F*A))*u;
cs_surf_neg = min(cs_surf_neg_ub, max(cs_surf_neg_lb, cs_surf_neg)); % saturation

cs_avg_pos = cs_max_p*(theta_min_pos + ((x(3) - theta_min_neg*cs_max_n)./((theta_max_neg - theta_min_neg)*cs_max_n))*(theta_max_pos - theta_min_pos));  % average concentration in pos electrodue (analytical expn using conservation of Li)
cs_surf_pos = cs_avg_pos + (8*R_pos/35)*x(1) + (R_pos/(35*Ds_pos*a_pos*L_pos*F*A))*u; % surface concentration of pos electrode
cs_surf_pos = min(cs_surf_pos_ub, max(cs_surf_pos_lb, cs_surf_pos)); % saturation

%% Compute overpotentials
eta_p = (2*R*T/F)*asinh(-u/(2*A*F*a_pos*L_pos*k_p*sqrt(ce*cs_surf_pos*(cs_max_p - cs_surf_pos))));
eta_n = (2*R*T/F)*asinh(u/(2*A*F*a_neg*L_neg*k_n*sqrt(ce*cs_surf_neg*(cs_max_n - cs_surf_neg))));

%% Compute OCPs
surf_theta_p  = cs_surf_pos/param.cs_max_p;
surf_theta_n  = cs_surf_neg/param.cs_max_n;

U_p = compute_Uocp_pos(surf_theta_p);
U_n = compute_Uocp_neg(surf_theta_n);

phi_pos = eta_p + U_p;
phi_neg = eta_n + U_n;

v_cell = phi_pos - phi_neg;

end
