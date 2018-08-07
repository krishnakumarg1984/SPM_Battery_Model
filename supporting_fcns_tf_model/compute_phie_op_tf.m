t_plus = ce_tf_subsys_params.tplus;
R = ce_tf_subsys_params.R;
T = ce_tf_subsys_params.Tref;
A = ce_tf_subsys_params.A;
F = ce_tf_subsys_params.F;
ce_init = ce_tf_subsys_params.ce_init;

phie_op_vector_results = nan(size(ce_neg_cc_tf));
phie_op_term1_vector = nan(size(ce_neg_cc_tf));
phie_op_term2_vector = nan(size(ce_neg_cc_tf));
 
progressbarText(0);
step_tf = 1;
while step_tf <= no_of_steps
    I_load_t = I_load_vector(step_tf);
    
    phie_op_term1 = (1-t_plus)*(2*R*T/F)*log(ce_pos_cc_tf(step_tf)./ce_neg_cc_tf(step_tf)); 
    phie_op_term1_vector(step_tf) = phie_op_term1;
%     kappa_eff_neg = electrolyteConductivity(ce_neg_cc_tf(step_tf), T, ce_tf_subsys_params, 'n');
%     kappa_eff_sep = electrolyteConductivity(ce_sep_tf_0p5Ls(step_tf), T, ce_tf_subsys_params, 's');
%     kappa_eff_pos = electrolyteConductivity(ce_pos_cc_tf(step_tf), T, ce_tf_subsys_params, 'p');
    
    kappa_eff_neg = electrolyteConductivity(ce_init, T, ce_tf_subsys_params, 'n');
    kappa_eff_sep = electrolyteConductivity(ce_init, T, ce_tf_subsys_params, 's');
    kappa_eff_pos = electrolyteConductivity(ce_init, T, ce_tf_subsys_params, 'p');

    phie_op_term2a = (-I_load_t/(2*A))*(Ln/kappa_eff_neg);
    phie_op_term2b = (-I_load_t/(2*A))*(2*(Ls/kappa_eff_sep));
    phie_op_term2c = (-I_load_t/(2*A))*(Lp/kappa_eff_pos);
    phie_op_term2 = phie_op_term2a + phie_op_term2b + phie_op_term2c;
	phie_op_term2_vector(step_tf) = phie_op_term2;
    
    phie_op_vector_results(step_tf) = phie_op_term1 + phie_op_term2;
    
    t0_local_tf = t0_local_tf + Ts;
    tf_local_tf = tf_local_tf + Ts;
    step_tf = step_tf + 1;
    progressbarText((step_tf-1)/no_of_steps);
end

clear t_plus R T A F kappa_eff_neg kappa_eff_sep kappa_eff_pos phie_op_term2 phie_op_term1 I_load_t;
