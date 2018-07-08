function ce_sep = compute_ce_sep_od_poly(a4,a3,a2,z,Ls,ce_init)
    ce_sep = a4.*((z/Ls).^2) + a3.*(z/Ls) + a2 + ce_init;
end