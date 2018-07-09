
Nn = 100;
Ns = Nn;
Np = Nn;

z_neg = linspace(0,Ln,Nn);
z_sep = linspace(0,Ls,Ns);
z_pos = linspace(0,Lp,Np);
%% Post-Process quadratic Results
a0_quadratic_vector = a_quadratic_vector_results(1,:);
a2_quadratic_vector = a_quadratic_vector_results(2,:);
a3_quadratic_vector = a_quadratic_vector_results(3,:);
a4_quadratic_vector = a_quadratic_vector_results(4,:);
a5_quadratic_vector = a_quadratic_vector_results(5,:);
a6_quadratic_vector = a_quadratic_vector_results(6,:);
a8_quadratic_vector = a_quadratic_vector_results(7,:);

% Assemble quadratic results at a given spatial point for all time
ce_neg_cc_quadratic = compute_ce_neg_quadratic(a2_quadratic_vector,a0_quadratic_vector,0);
ce_neg_sep_quadratic = compute_ce_neg_quadratic(a2_quadratic_vector,a0_quadratic_vector,Ln);

ce_pos_cc_quadratic = compute_ce_pos_quadratic(a8_quadratic_vector,a6_quadratic_vector,0);
ce_pos_sep_quadratic = compute_ce_pos_quadratic(a8_quadratic_vector,a6_quadratic_vector,Lp);

%%
run('ce_setup_plot_time_snapshots.m');

[~, tplot1_idx] = min(abs(time_vector_p2d-t_plot1));
if ~isempty(tplot1_idx)
    ce_neg_quadratic_tplot1 = compute_ce_neg_quadratic(a2_quadratic_vector(tplot1_idx),a0_quadratic_vector(tplot1_idx),z_neg);
    ce_sep_quadratic_tplot1 = compute_ce_sep_quadratic(a5_quadratic_vector(tplot1_idx),a4_quadratic_vector(tplot1_idx),a3_quadratic_vector(tplot1_idx),z_sep);
    ce_pos_quadratic_tplot1 = compute_ce_pos_quadratic(a8_quadratic_vector(tplot1_idx),a6_quadratic_vector(tplot1_idx),z_pos);
    ce_quadratic_tplot1 = [ce_neg_quadratic_tplot1,(ce_sep_quadratic_tplot1),fliplr(ce_pos_quadratic_tplot1)];
end

[~, tplot2_idx] = min(abs(time_vector_p2d-t_plot2));
if ~isempty(tplot2_idx)
    ce_neg_quadratic_tplot2 = compute_ce_neg_quadratic(a2_quadratic_vector(tplot2_idx),a0_quadratic_vector(tplot2_idx),z_neg);
    ce_sep_quadratic_tplot2 = compute_ce_sep_quadratic(a5_quadratic_vector(tplot2_idx),a4_quadratic_vector(tplot2_idx),a3_quadratic_vector(tplot2_idx),z_sep);
    ce_pos_quadratic_tplot2 = compute_ce_pos_quadratic(a8_quadratic_vector(tplot2_idx),a6_quadratic_vector(tplot2_idx),z_pos);
    ce_quadratic_tplot2 = [ce_neg_quadratic_tplot2,(ce_sep_quadratic_tplot2),fliplr(ce_pos_quadratic_tplot2)];
end

[~, tplot3_idx] = min(abs(time_vector_p2d-t_plot3));
if ~isempty(tplot3_idx)
    ce_neg_quadratic_tplot3 = compute_ce_neg_quadratic(a2_quadratic_vector(tplot3_idx),a0_quadratic_vector(tplot3_idx),z_neg);
    ce_sep_quadratic_tplot3 = compute_ce_sep_quadratic(a5_quadratic_vector(tplot3_idx),a4_quadratic_vector(tplot3_idx),a3_quadratic_vector(tplot3_idx),z_sep);
    ce_pos_quadratic_tplot3 = compute_ce_pos_quadratic(a8_quadratic_vector(tplot3_idx),a6_quadratic_vector(tplot3_idx),z_pos);
    ce_quadratic_tplot3 = [ce_neg_quadratic_tplot3,(ce_sep_quadratic_tplot3),fliplr(ce_pos_quadratic_tplot3)];
end

[~, tplot4_idx] = min(abs(time_vector_p2d-t_plot4));
if ~isempty(tplot4_idx)
    ce_neg_quadratic_tplot4 = compute_ce_neg_quadratic(a2_quadratic_vector(tplot4_idx),a0_quadratic_vector(tplot4_idx),z_neg);
    ce_sep_quadratic_tplot4 = compute_ce_sep_quadratic(a5_quadratic_vector(tplot4_idx),a4_quadratic_vector(tplot4_idx),a3_quadratic_vector(tplot4_idx),z_sep);
    ce_pos_quadratic_tplot4 = compute_ce_pos_quadratic(a8_quadratic_vector(tplot4_idx),a6_quadratic_vector(tplot4_idx),z_pos);
    ce_quadratic_tplot4 = [ce_neg_quadratic_tplot4,(ce_sep_quadratic_tplot4),fliplr(ce_pos_quadratic_tplot4)];
end

[~, tplot5_idx] = min(abs(time_vector_p2d-t_plot5));
if ~isempty(tplot5_idx)
    ce_neg_quadratic_tplot5 = compute_ce_neg_quadratic(a2_quadratic_vector(tplot5_idx),a0_quadratic_vector(tplot5_idx),z_neg);
    ce_sep_quadratic_tplot5 = compute_ce_sep_quadratic(a5_quadratic_vector(tplot5_idx),a4_quadratic_vector(tplot5_idx),a3_quadratic_vector(tplot5_idx),z_sep);
    ce_pos_quadratic_tplot5 = compute_ce_pos_quadratic(a8_quadratic_vector(tplot5_idx),a6_quadratic_vector(tplot5_idx),z_pos);
    ce_quadratic_tplot5 = [ce_neg_quadratic_tplot5,(ce_sep_quadratic_tplot5),fliplr(ce_pos_quadratic_tplot5)];
end

[~, tplot6_idx] = min(abs(time_vector_p2d-t_plot6));
if ~isempty(tplot6_idx)
    ce_neg_quadratic_tplot6 = compute_ce_neg_quadratic(a2_quadratic_vector(tplot6_idx),a0_quadratic_vector(tplot6_idx),z_neg);
    ce_sep_quadratic_tplot6 = compute_ce_sep_quadratic(a5_quadratic_vector(tplot6_idx),a4_quadratic_vector(tplot6_idx),a3_quadratic_vector(tplot6_idx),z_sep);
    ce_pos_quadratic_tplot6 = compute_ce_pos_quadratic(a8_quadratic_vector(tplot6_idx),a6_quadratic_vector(tplot6_idx),z_pos);
    ce_quadratic_tplot6 = [ce_neg_quadratic_tplot6,(ce_sep_quadratic_tplot6),fliplr(ce_pos_quadratic_tplot6)];
end
