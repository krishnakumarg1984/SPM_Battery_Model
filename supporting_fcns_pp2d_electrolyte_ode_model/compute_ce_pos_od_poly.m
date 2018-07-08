function ce_pos_poly = compute_ce_pos_od_poly(a6,a5,z,Lp,ce_init)
    ce_pos_poly = a6.*((z./Lp).^2) + a5 + ce_init;
end