function mean_free_path = get_v_dependent_mean_free_path(v, sigma)
%GET_V_DEPENDENT_MEAN_FREE_PATH 

bulk_density_helium = 0.0219; % Angström ^-3
density_droplet = 0.8*bulk_density_helium;




global sigma_dependent_on_v

if sigma_dependent_on_v
sigma_v = sigma*0.6* exp(- v.^2 / 8^2) + 0.4*sigma;


mean_free_path= 1./(density_droplet*  sigma_v);
else

    mean_free_path = 1./(density_droplet*  sigma);
end


% debug
% debug
 % v = 0.1:0.1:20
 % hold on
 % sigma_v = sigma*0.8* exp(- v.^2 / 12^2) + 0.2*sigma;
 % plot(v, sigma_v);

% three new parameters here:
% sigma max, sigma min, decay constant (slope of log(sigma(E)) over log(E))

end

