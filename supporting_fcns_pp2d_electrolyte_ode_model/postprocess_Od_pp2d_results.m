z_neg = linspace(0,Ln,extradata_for_solver.pp2d_model_params.Nn);
z_sep = linspace(0,Ls,extradata_for_solver.pp2d_model_params.Ns);
z_pos = linspace(0,Lp,extradata_for_solver.pp2d_model_params.Np);

%% Post-Process Polynomial Results

a0_od_poly_solved_alltime = a_od_poly_vector_results(1,:);
a1_od_poly_solved_alltime = a_od_poly_vector_results(2,:);
a2_od_poly_solved_alltime = a_od_poly_vector_results(3,:);
a3_od_poly_solved_alltime = a_od_poly_vector_results(4,:);
a4_od_poly_solved_alltime = a_od_poly_vector_results(5,:);
a5_od_poly_solved_alltime = a_od_poly_vector_results(6,:);
a6_od_poly_solved_alltime = a_od_poly_vector_results(7,:);

ce_neg_cc_od_poly = compute_ce_neg_od_poly(a1_od_poly_solved_alltime,a0_od_poly_solved_alltime,0,Ln,ce_init);
