mu = @(N) -7.21 + 17.71*N.^(-1/3) - 5.95*N.^(-2/3); 

T = 300; 
N = 1000:100:20000;

E = (T*3/2)*k_B;

D0 = 1.556*eV;
E_pump = 1.9*eV;

E_pump = 1.9*eV*2;
D1 = 2.526*eV;

E_lim =2*127*u/2*(40)^2;


E = E_pump - D1 - E_lim;

num_evap =  E  ./(mu(N)*k_B);


plot(N, num_evap./N*100)
plot(N, num_evap)