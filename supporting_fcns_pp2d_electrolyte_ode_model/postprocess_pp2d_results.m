z_neg = linspace(0,Ln,extradata_for_solver.pp2d_model_params.Nn);
z_sep = linspace(0,Ls,extradata_for_solver.pp2d_model_params.Ns);
z_pos = linspace(0,Lp,extradata_for_solver.pp2d_model_params.Np);

%% Post-Process Polynomial Results

a0_poly_solved_alltime = a_poly_vector_results(1,:);
a1_poly_solved_alltime = a_poly_vector_results(2,:);
a2_poly_solved_alltime = a_poly_vector_results(3,:);
a3_poly_solved_alltime = a_poly_vector_results(4,:);
a4_poly_solved_alltime = a_poly_vector_results(5,:);
a5_poly_solved_alltime = a_poly_vector_results(6,:);
a6_poly_solved_alltime = a_poly_vector_results(7,:);

% Assemble P2D results at a given spatial point for all time
ce_neg_cc_poly = compute_ce_neg_poly(a1_poly_solved_alltime,a0_poly_solved_alltime,0);
ce_neg_sep_poly = compute_ce_neg_poly(a1_poly_solved_alltime,a0_poly_solved_alltime,Ln);

ce_pos_cc_poly = compute_ce_pos_poly(a6_poly_solved_alltime,a5_poly_solved_alltime,0);
ce_pos_sep_poly = compute_ce_pos_poly(a6_poly_solved_alltime,a5_poly_solved_alltime,Lp);

[~, pp2d_tplot1_idx] = min(abs(time_p2d-t_plot1));
if ~isempty(pp2d_tplot1_idx)
    ce_neg_poly_tplot1 = compute_ce_neg_poly(a1_poly_solved_alltime(pp2d_tplot1_idx),a0_poly_solved_alltime(pp2d_tplot1_idx),z_neg);
    ce_sep_poly_tplot1 = compute_ce_sep_poly(a4_poly_solved_alltime(pp2d_tplot1_idx),a3_poly_solved_alltime(pp2d_tplot1_idx),a2_poly_solved_alltime(pp2d_tplot1_idx),z_sep);
    ce_pos_poly_tplot1 = compute_ce_pos_poly(a6_poly_solved_alltime(pp2d_tplot1_idx),a5_poly_solved_alltime(pp2d_tplot1_idx),z_pos);
    ce_poly_tplot1 = [ce_pos_poly_tplot1,fliplr(ce_sep_poly_tplot1),fliplr(ce_neg_poly_tplot1)];
end

[~, pp2d_tplot2_idx] = min(abs(time_p2d-t_plot2));
if ~isempty(pp2d_tplot2_idx)
    ce_neg_poly_tplot2 = compute_ce_neg_poly(a1_poly_solved_alltime(pp2d_tplot2_idx),a0_poly_solved_alltime(pp2d_tplot2_idx),z_neg);
    ce_sep_poly_tplot2 = compute_ce_sep_poly(a4_poly_solved_alltime(pp2d_tplot2_idx),a3_poly_solved_alltime(pp2d_tplot2_idx),a2_poly_solved_alltime(pp2d_tplot2_idx),z_sep);
    ce_pos_poly_tplot2 = compute_ce_pos_poly(a6_poly_solved_alltime(pp2d_tplot2_idx),a5_poly_solved_alltime(pp2d_tplot2_idx),z_pos);
    ce_poly_tplot2 = [ce_pos_poly_tplot2,fliplr(ce_sep_poly_tplot2),fliplr(ce_neg_poly_tplot2)];
end

[~, pp2d_tplot3_idx] = min(abs(time_p2d-t_plot3));
if ~isempty(pp2d_tplot3_idx)
    ce_neg_poly_tplot3 = compute_ce_neg_poly(a1_poly_solved_alltime(pp2d_tplot3_idx),a0_poly_solved_alltime(pp2d_tplot3_idx),z_neg);
    ce_sep_poly_tplot3 = compute_ce_sep_poly(a4_poly_solved_alltime(pp2d_tplot3_idx),a3_poly_solved_alltime(pp2d_tplot3_idx),a2_poly_solved_alltime(pp2d_tplot3_idx),z_sep);
    ce_pos_poly_tplot3 = compute_ce_pos_poly(a6_poly_solved_alltime(pp2d_tplot3_idx),a5_poly_solved_alltime(pp2d_tplot3_idx),z_pos);
    ce_poly_tplot3 = [ce_pos_poly_tplot3,fliplr(ce_sep_poly_tplot3),fliplr(ce_neg_poly_tplot3)];
end

[~, pp2d_tplot4_idx] = min(abs(time_p2d-t_plot4));
if ~isempty(pp2d_tplot4_idx)
    ce_neg_poly_tplot4 = compute_ce_neg_poly(a1_poly_solved_alltime(pp2d_tplot4_idx),a0_poly_solved_alltime(pp2d_tplot4_idx),z_neg);
    ce_sep_poly_tplot4 = compute_ce_sep_poly(a4_poly_solved_alltime(pp2d_tplot4_idx),a3_poly_solved_alltime(pp2d_tplot4_idx),a2_poly_solved_alltime(pp2d_tplot4_idx),z_sep);
    ce_pos_poly_tplot4 = compute_ce_pos_poly(a6_poly_solved_alltime(pp2d_tplot4_idx),a5_poly_solved_alltime(pp2d_tplot4_idx),z_pos);
    ce_poly_tplot4 = [ce_pos_poly_tplot4,fliplr(ce_sep_poly_tplot4),fliplr(ce_neg_poly_tplot4)];
end

[~, pp2d_tplot5_idx] = min(abs(time_p2d-t_plot5));
if ~isempty(pp2d_tplot5_idx)
    ce_neg_poly_tplot5 = compute_ce_neg_poly(a1_poly_solved_alltime(pp2d_tplot5_idx),a0_poly_solved_alltime(pp2d_tplot5_idx),z_neg);
    ce_sep_poly_tplot5 = compute_ce_sep_poly(a4_poly_solved_alltime(pp2d_tplot5_idx),a3_poly_solved_alltime(pp2d_tplot5_idx),a2_poly_solved_alltime(pp2d_tplot5_idx),z_sep);
    ce_pos_poly_tplot5 = compute_ce_pos_poly(a6_poly_solved_alltime(pp2d_tplot5_idx),a5_poly_solved_alltime(pp2d_tplot5_idx),z_pos);
    ce_poly_tplot5 = [ce_pos_poly_tplot5,fliplr(ce_sep_poly_tplot5),fliplr(ce_neg_poly_tplot5)];
end

[~, pp2d_tplot6_idx] = min(abs(time_p2d-t_plot6));
if ~isempty(pp2d_tplot6_idx)
    ce_neg_poly_tplot6 = compute_ce_neg_poly(a1_poly_solved_alltime(pp2d_tplot6_idx),a0_poly_solved_alltime(pp2d_tplot6_idx),z_neg);
    ce_sep_poly_tplot6 = compute_ce_sep_poly(a4_poly_solved_alltime(pp2d_tplot6_idx),a3_poly_solved_alltime(pp2d_tplot6_idx),a2_poly_solved_alltime(pp2d_tplot6_idx),z_sep);
    ce_pos_poly_tplot6 = compute_ce_pos_poly(a6_poly_solved_alltime(pp2d_tplot6_idx),a5_poly_solved_alltime(pp2d_tplot6_idx),z_pos);
    ce_poly_tplot6 = [ce_pos_poly_tplot6,fliplr(ce_sep_poly_tplot6),fliplr(ce_neg_poly_tplot6)];
end
