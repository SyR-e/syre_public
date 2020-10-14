%
fdfq_idiq [T, Fd, Fq, L, pf] = ...
         compute_wk_point2_(id, iq, i, p_Ld0, p_G0Fp, p_G0F, p_Lq0, p_F0Gp, p_F0G, F0Gend, G0Fend, T_cost)
%
Fd = polyval(p_Ld0,id)-(polyval(p_G0Fp,id) .*  (polyval(p_F0G , iq))) / F0Gend;
Fq = polyval(p_Lq0,iq)-(polyval(p_F0Gp, iq).*  (polyval(p_G0F , id))) / G0Fend;
L  = abs(Fd + j *Fq);
%
T = Fd .* iq - Fq .* id;
pf = T ./ ( L .* i);
T = T_cost * T;
