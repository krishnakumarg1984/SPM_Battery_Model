function ce_sep = compute_ce_sep_tf(a5,a4,a3,z)
    ce_sep = a5.*(z.^2) + a4*z + a3;
end