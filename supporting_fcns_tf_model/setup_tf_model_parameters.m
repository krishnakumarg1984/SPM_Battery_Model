%% Load the identified Qe transfer functions
load('Qen_tf_4p3z_scaled.mat');
load('Qes_tf_2p1z_scaled.mat');
load('Qep_tf_4p3z_scaled.mat');

tf_temp = tf; clear tf;
Qen_tf = tf(Qen_tfest_scaled);
Qes_tf = tf(Qes_tfest_scaled);
Qep_tf = tf(Qep_tfest_scaled);
tf = tf_temp; clear tf_temp;

tf_scale_factor = 1e3;

ce_tf_subsys_params = define_ce_tf_model_parameters(soc_init_pct);
eps_n = ce_tf_subsys_params.eps_n;
eps_s = ce_tf_subsys_params.eps_s;
eps_p = ce_tf_subsys_params.eps_p;
ce_init = ce_tf_subsys_params.ce_init;
Ln = ce_tf_subsys_params.len_n;
Ls = ce_tf_subsys_params.len_s;
Lp = ce_tf_subsys_params.len_p;

Qe_ML_vector_init = [eps_n*ce_init*Ln;eps_s*ce_init*Ls;eps_p*ce_init*Lp]; % at t = 0

% Qe_ML_vector_results = nan(length(Qe_ML_vector_init),no_of_steps);
% Qe_ML_vector_results(:,1) = Qe_ML_vector_init;
no_of_steps = ceil(t_finish/Ts) + 1; % max no. of steps (assuming no cutoff)
tf_interp_vector = 0:Ts:t_finish;
I_load_vector = I_1C*(interp1(C_rate_profile(:,1),C_rate_profile(:,2),tf_interp_vector,'previous','extrap'));

del_Qen_ML = (1/tf_scale_factor)*lsim(Qen_tf,I_load_vector,tf_interp_vector,[],'zoh');
del_Qes_ML = (1/tf_scale_factor)*lsim(Qes_tf,I_load_vector,tf_interp_vector,[],'zoh');
del_Qep_ML = (1/tf_scale_factor)*lsim(Qep_tf,I_load_vector,tf_interp_vector,[],'zoh');

Qe_tf_vector_results = [del_Qen_ML';del_Qes_ML';del_Qep_ML'] + Qe_ML_vector_init;