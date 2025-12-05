function [N, R] = get_dropletsize(p, T,d)
%GET_DROPLETSIZE calculate the mean droplet size of helium droplets from
%source parameters
% 
% if size(p,2)> size(p,1)
%     p = p';
% end
% if size(T, 2)<size(T,1)
%     T= T';
% end



if nargin<3
    d = 5; % nozzle diameter in micrometer
end

k1 = 4E5;
k2 = 0.97;
k3 = -3.88;
k4 = 2;


N = k1*(p.^k2).*(T.^k3).*d.^k4;

rho_bulk = 0.0218;
rho_droplet = rho_bulk*0.8;
R = (3/(4*pi*rho_droplet))^(1/3) * N.^(1/3);

%R = 2.22*N.^(1/3); %radius in angström

fprintf('R = %.2f Angström ---- ', R);

fprintf('sigma = %.2f Gb \n', mean(pi*R.^2)*0.1);% crossection in gigabarn

end

