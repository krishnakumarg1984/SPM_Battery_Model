t_plus = ce_tf_subsys_params.tplus;
R = ce_tf_subsys_params.R;
T = ce_tf_subsys_params.Tref;
A = ce_tf_subsys_params.A;
F = ce_tf_subsys_params.F;
% return;

progressbarText(0);
step_tf = 1;
while step_tf <= no_of_steps
    I_load_t = I_load_vector(step_tf);
    
    phie_op_term1 = (1-t_plus)*(2*R*T/F)*log(ce_pos_cc_tf(step_tf)./ce_neg_cc_tf(step_tf)); 
    
    kappa_eff_neg = electrolyteConductivity(ce_pos_cc_tf(step_tf), T, ce_tf_subsys_params, 'n');

    t0_local_tf = t0_local_tf + Ts;
    tf_local_tf = tf_local_tf + Ts;
    step_tf = step_tf + 1;
    progressbarText((step_tf-1)/no_of_steps);
end


ce_tf_subsys_params.De_eff_pos = electrolyteDiffusionCoefficients(params_temp.ce_init,params_temp.Tref,params_temp,'p');


phie_opt_term2 = I_load_vector./(2*A)
clear t_plus R T;
