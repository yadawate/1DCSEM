function Ratm = loseth_tm(krho)
global k h srln n_int zz eps_h Srdpt
Rm_tm = 0;
for kd = n_int:-1:srln+1
     rtm_n = (eps_h(kd) .* kzf(krho,k(kd+1))) - (eps_h(kd+1) .* kzf(krho,k(kd)));
     rtm_d = (eps_h(kd) .* kzf(krho,k(kd+1))) + (eps_h(kd+1) .* kzf(krho,k(kd)));
     rtm = rtm_n ./ rtm_d;
     Rm_tm = (Rm_tm + rtm) .* (exp(2j .* kzf(krho,k(kd)) .* h(kd))) ./ (1 + (rtm .* Rm_tm));
end
 
rtm_n = (eps_h(2) .* kzf(krho,k(3))) - (eps_h(3) .* kzf(krho,k(2)));
rtm_d = (eps_h(2) .* kzf(krho,k(3))) + (eps_h(3) .* kzf(krho,k(2)));
rtm = rtm_n ./ rtm_d;
Rm_tm = Rm_tm + rtm.* exp(1j .* kzf(krho,k(srln)) .* (2*(h(srln)-Srdpt)));

rtm_n = eps_h(1) .* kzf(krho,k(2)) - eps_h(2) .* kzf(krho,k(1));
rtm_d = eps_h(1) .* kzf(krho,k(2)) + eps_h(2) .* kzf(krho,k(1));
rtm = rtm_n ./ rtm_d;
Rm_tmu = -rtm.*exp(1j.*kzf(krho,k(2)).*(2*Srdpt)); 

n1m = Rm_tm .* (1 + Rm_tmu) .* exp(1j .* kzf(krho,k(srln)) .* zz);
n2m = Rm_tmu .* (1 + Rm_tm) .* exp(1j .* kzf(krho,k(srln)) .* zz);
dm = 1 - (Rm_tm .* Rm_tmu);
Ratm = (n1m + n2m) ./ dm;

end
