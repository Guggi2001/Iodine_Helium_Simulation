function r_accepted_all_radii = generate_radial_samples_3d(droplet_radii, droplet_potential,steepness, binding_energy, E_max, T_droplet, r_step, debug_plot)
%GENERATE_RADIAL_SAMPLES generate random samples of molecules inside helium
%droplet, given a droplet potential function and temperature
eV = 1.602e-19; %joule
boltzmann_constant = 0.08617*1E-3*eV; % joule per kelvin

p_boltzmann = @(EE, TT, ZZ) exp(-EE./(boltzmann_constant*TT))/ZZ;
E = [0:0.001:E_max]/1000*eV; % energies in joule

normalization = trapz(E, p_boltzmann(E, T_droplet, 1)); 

p = @(E) p_boltzmann(E, T_droplet, normalization);


if size(droplet_radii,2)>1 % flip if input array has wrong orientation

    droplet_radii = droplet_radii';
end

unique_radii = unique(droplet_radii);

r_accepted_all_radii = [];

for droplet_radius =  unique_radii'
    num_samples = sum(droplet_radius==droplet_radii);

    r_max = droplet_radius*2;

r = 0:r_step:r_max;
normalization = trapz( r,  p(droplet_potential([steepness, binding_energy, droplet_radius], r)/1000*eV).*r.^2); % this needs to include sin(theta) i think

p_radius = @(r) p(droplet_potential([steepness, binding_energy, droplet_radius],r)/1000*eV).*r.^2/normalization;



y_max = max(p_radius(r));


r_accepted_total = [];

while length(r_accepted_total) < num_samples
% random sampling radii from this distribution
num_proposals = 1000;


r_proposal = rand(num_proposals, 1)*r_max;

y = p_radius(r_proposal);

trial_probability = rand(num_proposals,1)*y_max;
accept = trial_probability<y;

r_accepted = r_proposal(accept);
r_accepted_total = [r_accepted_total; r_accepted];

accept_rate = sum(accept)/size(r_proposal,1);

end




if debug_plot
    if length(unique_radii)<5
% 
% figure
% plot(E/eV*1000, p(E));
% xlim([0,30]);

figure
plot(r, p_radius(r));
hold on
bayes_hist(r_accepted_total, [0,r_max], false, 'red');
title(sprintf('acceptance rate: %.2f\n', accept_rate));
xlabel('radius / Anström')
ylabel('probabilitiy');

r_test = [0:r_step:100];
pot = droplet_potential([steepness, binding_energy, droplet_radius],r_test);
yyaxis right
plot(r_test, pot/1000*eV/boltzmann_constant);
hline(0.34);
ylabel('T / K');
legend('p(r)', 'histogram of sampled radii', 'droplet binding energy')
%set(gca,'YScale', 'log');
set(gca,'Fontsize', 15);

    else
       
    end

end


r_accepted_all_radii = [r_accepted_all_radii; r_accepted_total(1:num_samples)];



if debug_plot
%bayes_hist(r_accepted_total, [0,r_max], false, 'red');
end

end


end

