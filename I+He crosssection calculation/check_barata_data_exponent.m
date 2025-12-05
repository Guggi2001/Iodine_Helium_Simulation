load("barata_data_XeHe.mat");

plot(E, Qs);

fun = @(beta, E) beta(1)*E.^(beta(2));

beta = nlinfit(E, Qs, fun, [1,-1]);
hold on
plot(E, fun(beta,E));


E_kin = E*eV/2;

v = sqrt(2*E_kin/(40*u));

figure
 plot(v, fun(beta,E));
 