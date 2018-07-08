function [res_XZ, flag, new_data] = pp2d_electrolyte_model_residuals(t,XZ,dXZ_dt,ida_user_data)

flag     = 0;  % These two flags are not used but required by IDA(s) solver
new_data = []; % These two flags are not used but required by IDA(s) solver

t0 = ida_user_data.t0;
tf = ida_user_data.tf;

ida_user_data.I = ida_user_data.load_curr_anonymousFcn(t,t0,tf); % the RHS evaluates the applied load current at present solver-time 't'

res_odes_pp2d_electrolyte = ode_res_pp2d_electrolyte_model(t,XZ,dXZ_dt,ida_user_data); % residuals of ODE eqns
res_alg_pp2d_electrolyte = alg_res_pp2d_electrolyte_model(XZ,ida_user_data); % residuals of algebraic eqns
res_XZ = [res_odes_pp2d_electrolyte;res_alg_pp2d_electrolyte]; % concatenate to form overall residual vector

end
