function ce_neg_poly = compute_ce_neg_od_poly(a1,a0,z,Ln,ce_init)
    ce_neg_poly = a1.*((z./Ln).^2) + a0 + ce_init;
end