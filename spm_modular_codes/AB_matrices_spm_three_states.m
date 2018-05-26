function [A_cts, A_disc, B_cts, B_disc] = AB_matrices_spm_three_states(param,Ts)
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

%% System & Input matrices for continuous and discrete-time implementations

A_cts = [-30*Ds_pos/(R_pos^2),                    0, 0; ...
                            0, -30*Ds_neg/(R_neg^2), 0; ...
                            0,                    0, 0];

A_disc = expm(A_cts*Ts);

B_cts = [ (45/2)/(R_pos^2*a_pos*L_pos*F*A); ...
         (-45/2)/(R_neg^2*a_neg*L_neg*F*A); ...
             (-3/(R_neg*a_neg*L_neg*F*A))];

B_disc = nan(size(B_cts));
B_disc(1) = B_cts(1)*(exp(A_cts(1,1)*Ts)-1)/A_cts(1,1);
B_disc(2) = B_cts(2)*(exp(A_cts(2,2)*Ts)-1)/A_cts(2,2);
B_disc(3) = B_cts(3)*Ts;

% Through built-in 'c2d' command by the dummy system-method
% (tricking MATLAB to believe we have an LTI system).
% The output matrix,C is chosen so that the states-themselves are the
% outputs with no feed-through term involved.

C_dummy = [1, 1, 1]; D_dummy = 0;
sys_cts = ss(A_cts,B_cts,C_dummy,D_dummy);
sys_disc = c2d(sys_cts,Ts);

% A_disc = sys_disc.A
% B_disc = sys_disc.B

end
