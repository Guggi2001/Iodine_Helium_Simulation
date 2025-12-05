function y = droplet_potential(beta,x)

y = ((erf((x-beta(3))/beta(1))*1+1)/2)*beta(2);


end

