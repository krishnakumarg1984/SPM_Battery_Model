
Nn = 100;
Ns = Nn;
Np = Nn;

z_neg = linspace(0,Ln,Nn);
z_sep = linspace(0,Ls,Ns);
z_pos = linspace(0,Lp,Np);

%% Post-Process tf Results
a0_tf_vector = a_tf_vector_results(1,:);
a2_tf_vector = a_tf_vector_results(2,:);
a3_tf_vector = a_tf_vector_results(3,:);
a4_tf_vector = a_tf_vector_results(4,:);
a5_tf_vector = a_tf_vector_results(5,:);
a6_tf_vector = a_tf_vector_results(6,:);
a8_tf_vector = a_tf_vector_results(7,:);

% Assemble tf results at a given spatial point for all time
ce_neg_cc_tf = compute_ce_neg_tf(a2_tf_vector,a0_tf_vector,0);
ce_neg_sep_tf = compute_ce_neg_tf(a2_tf_vector,a0_tf_vector,Ln);

% ce_sep_tf_0 = compute_ce_sep_tf(a5_tf_vector,a4_tf_vector,a3_tf_vector,0);
ce_sep_tf_0p5Ls = compute_ce_sep_tf(a5_tf_vector,a4_tf_vector,a3_tf_vector,0.5*Ls);
% ce_sep_tf_Ls = compute_ce_sep_tf(a5_tf_vector,a4_tf_vector,a3_tf_vector,Ls);

ce_pos_cc_tf = compute_ce_pos_tf(a8_tf_vector,a6_tf_vector,0);
ce_pos_sep_tf = compute_ce_pos_tf(a8_tf_vector,a6_tf_vector,Lp);

%%
run('ce_setup_plot_time_snapshots.m');

[~, tplot1_idx] = min(abs(time_vector_p2d-t_plot1));
if ~isempty(tplot1_idx)
    ce_neg_tf_tplot1 = compute_ce_neg_tf(a2_tf_vector(tplot1_idx),a0_tf_vector(tplot1_idx),z_neg);
    ce_sep_tf_tplot1 = compute_ce_sep_tf(a5_tf_vector(tplot1_idx),a4_tf_vector(tplot1_idx),a3_tf_vector(tplot1_idx),z_sep);
    ce_pos_tf_tplot1 = compute_ce_pos_tf(a8_tf_vector(tplot1_idx),a6_tf_vector(tplot1_idx),z_pos);
    ce_tf_tplot1 = [ce_neg_tf_tplot1,(ce_sep_tf_tplot1),fliplr(ce_pos_tf_tplot1)];
end

[~, tplot2_idx] = min(abs(time_vector_p2d-t_plot2));
if ~isempty(tplot2_idx)
    ce_neg_tf_tplot2 = compute_ce_neg_tf(a2_tf_vector(tplot2_idx),a0_tf_vector(tplot2_idx),z_neg);
    ce_sep_tf_tplot2 = compute_ce_sep_tf(a5_tf_vector(tplot2_idx),a4_tf_vector(tplot2_idx),a3_tf_vector(tplot2_idx),z_sep);
    ce_pos_tf_tplot2 = compute_ce_pos_tf(a8_tf_vector(tplot2_idx),a6_tf_vector(tplot2_idx),z_pos);
    ce_tf_tplot2 = [ce_neg_tf_tplot2,(ce_sep_tf_tplot2),fliplr(ce_pos_tf_tplot2)];
end

[~, tplot3_idx] = min(abs(time_vector_p2d-t_plot3));
if ~isempty(tplot3_idx)
    ce_neg_tf_tplot3 = compute_ce_neg_tf(a2_tf_vector(tplot3_idx),a0_tf_vector(tplot3_idx),z_neg);
    ce_sep_tf_tplot3 = compute_ce_sep_tf(a5_tf_vector(tplot3_idx),a4_tf_vector(tplot3_idx),a3_tf_vector(tplot3_idx),z_sep);
    ce_pos_tf_tplot3 = compute_ce_pos_tf(a8_tf_vector(tplot3_idx),a6_tf_vector(tplot3_idx),z_pos);
    ce_tf_tplot3 = [ce_neg_tf_tplot3,(ce_sep_tf_tplot3),fliplr(ce_pos_tf_tplot3)];
end

[~, tplot4_idx] = min(abs(time_vector_p2d-t_plot4));
if ~isempty(tplot4_idx)
    ce_neg_tf_tplot4 = compute_ce_neg_tf(a2_tf_vector(tplot4_idx),a0_tf_vector(tplot4_idx),z_neg);
    ce_sep_tf_tplot4 = compute_ce_sep_tf(a5_tf_vector(tplot4_idx),a4_tf_vector(tplot4_idx),a3_tf_vector(tplot4_idx),z_sep);
    ce_pos_tf_tplot4 = compute_ce_pos_tf(a8_tf_vector(tplot4_idx),a6_tf_vector(tplot4_idx),z_pos);
    ce_tf_tplot4 = [ce_neg_tf_tplot4,(ce_sep_tf_tplot4),fliplr(ce_pos_tf_tplot4)];
end

[~, tplot5_idx] = min(abs(time_vector_p2d-t_plot5));
if ~isempty(tplot5_idx)
    ce_neg_tf_tplot5 = compute_ce_neg_tf(a2_tf_vector(tplot5_idx),a0_tf_vector(tplot5_idx),z_neg);
    ce_sep_tf_tplot5 = compute_ce_sep_tf(a5_tf_vector(tplot5_idx),a4_tf_vector(tplot5_idx),a3_tf_vector(tplot5_idx),z_sep);
    ce_pos_tf_tplot5 = compute_ce_pos_tf(a8_tf_vector(tplot5_idx),a6_tf_vector(tplot5_idx),z_pos);
    ce_tf_tplot5 = [ce_neg_tf_tplot5,(ce_sep_tf_tplot5),fliplr(ce_pos_tf_tplot5)];
end

[~, tplot6_idx] = min(abs(time_vector_p2d-t_plot6));
if ~isempty(tplot6_idx)
    ce_neg_tf_tplot6 = compute_ce_neg_tf(a2_tf_vector(tplot6_idx),a0_tf_vector(tplot6_idx),z_neg);
    ce_sep_tf_tplot6 = compute_ce_sep_tf(a5_tf_vector(tplot6_idx),a4_tf_vector(tplot6_idx),a3_tf_vector(tplot6_idx),z_sep);
    ce_pos_tf_tplot6 = compute_ce_pos_tf(a8_tf_vector(tplot6_idx),a6_tf_vector(tplot6_idx),z_pos);
    ce_tf_tplot6 = [ce_neg_tf_tplot6,(ce_sep_tf_tplot6),fliplr(ce_pos_tf_tplot6)];
end
