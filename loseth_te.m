function Rate = loseth_te(krho)
global k h srln n_int mu_h zz Srdpt
Rm_te = 0;
for kd = n_int:-1:srln+1
    rte_n = (mu_h(kd+1) .* kzf(krho,k(kd))) - (mu_h(kd) .* kzf(krho,k(kd+1)));
    rte_d = (mu_h(kd+1) .* kzf(krho,k(kd))) + (mu_h(kd) .* kzf(krho,k(kd+1)));
    rte = rte_n ./ rte_d;
    Rm_te = (Rm_te + rte) .* (exp(2j .* kzf(krho,k(kd)) .* h(kd)))./ (1 + (rte .* Rm_te));
end

rte_n = (mu_h(3) .* kzf(krho,k(2))) - (mu_h(2) .* kzf(krho,k(3)));
rte_d = (mu_h(3) .* kzf(krho,k(2))) + (mu_h(2) .* kzf(krho,k(3)));
rte = rte_n ./ rte_d;
Rm_te = Rm_te + rte.* exp(1j .* kzf(krho,k(srln)) .* (2*(h(srln)-Srdpt)));

rte_n = (mu_h(2) .* kzf(krho,k(1))) - (mu_h(1) .* kzf(krho,k(2)));
rte_d = (mu_h(2) .* kzf(krho,k(1))) + (mu_h(1) .* kzf(krho,k(2)));
rte = rte_n ./ rte_d;
Rm_teu = -rte.*exp(1j.*kzf(krho,k(2)).*(2*Srdpt));

n1e = Rm_te .* (1 + Rm_teu) .* exp(1j .* kzf(krho,k(srln)) .* zz);
n2e = Rm_teu .* (1 + Rm_te) .* exp(1j .* kzf(krho,k(srln)) .* zz);
de = 1 - (Rm_te .* Rm_teu);
Rate = (n1e + n2e) ./ de;
 
end