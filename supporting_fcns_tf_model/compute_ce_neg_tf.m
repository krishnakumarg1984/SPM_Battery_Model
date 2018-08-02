function ce_neg = compute_ce_neg_tf(a2,a0,z)
    ce_neg = a2.*(z.^2) + a0;
end