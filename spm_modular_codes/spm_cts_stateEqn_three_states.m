function dx_dt = spm_cts_stateEqn_three_states(x,u,param)
% returns dx_dt given the vector x and input u at any given time-step
% x(1) = q_pos, x(2) = q_neg, x(3) = cs_avg_neg
% u = load current (A). Positive implies discharge.

% Copyright (c) 2018 Gopalakrishnan, Krishnakumar <krishnak@vt.edu>
% Author: Gopalakrishnan, Krishnakumar <krishnak@vt.edu>

F = param.F;
A = param.A;

R_pos  = param.R_p;
R_neg  = param.R_n;

Ds_pos = param.D_p;
Ds_neg = param.D_n;

a_pos  = param.a_p;
a_neg  = param.a_n;

L_pos  = param.len_p;
L_neg  = param.len_n;

a_x1_cts = -30*Ds_pos/(R_pos^2);
b_x1_cts = (45/2)/(R_pos^2*a_pos*L_pos*F*A);

a_x2_cts = -30*Ds_neg/(R_neg^2);
b_x2_cts = (-45/2)/(R_neg^2*a_neg*L_neg*F*A);

dx1_dt = a_x1_cts*x(1) + b_x1_cts*u;
dx2_dt = a_x2_cts*x(2) + b_x2_cts*u;
dx3_dt = (-3/(R_neg*a_neg*L_neg*F*A))*u;

dx_dt  = [dx1_dt; dx2_dt ; dx3_dt];

end
