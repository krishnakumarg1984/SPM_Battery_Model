function ce_pos = compute_ce_pos_tf(a8,a6,z)
    ce_pos = a8.*z.^2 + a6;
end