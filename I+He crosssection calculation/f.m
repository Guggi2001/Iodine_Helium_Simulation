function [res, diagnostic] = f(v, theta, lmax, V)


% calculate crosssection according to 10.1109/NSSMIC.2004.1462242
rmin  =1;
rmax = 300;

mu = 127*4/(127 + 4); % in u

u = 1.66053907e-27; %kg
hbar_SI = 1.05457182E-34;
eV = 1.602e-19; %joule;
fs = 1E-15;


hbar = 0.6582845318; % eV femtosecond

unitfactor1 = (2*3.8*u*1*eV/(hbar_SI^2) / (1/1E-10)^2 )  /  (2*3.8*1/hbar^2);% convert u/(fs^2 eV) into 1/Angström^2

%https://www.wolframalpha.com/input?i=u+%2F%28eV+fs%5E2%29+%2F+%28+1%2FAngstr%C3%B6m%5E2%29
% proof of unitfactor:
% direct: https://www.wolframalpha.com/input?i=2*3.8+u+*+1+eV+%2F+%28hbar%5E2+%29+%2F+%281%2FAngstr%C3%B6m%5E2%29
% using unitfactor1: https://www.wolframalpha.com/input?i=2*3.8+*+1+%2F+%280.6582%5E2%29+*+103.642


res = zeros(length(v),length(theta));


unitfactor2 = ((3.8*u*1000/hbar_SI) /( 1/1E-10) )  / ( 3.8*10/hbar);% 

k_array =  mu *v/hbar * unitfactor2; % in 1/Angström



% determine r0 for each k (closest approach)
% rsample = linspace(1,30, 1000);
% for i=1:length(k)
%     K = k(i);
%     lsample = 6;
%     plot(rsample, K^2 - 2*mu*V(rsample)/hbar^2*unitfactor1- (lsample+1/2)^2./rsample.^2);
%     yyaxis right
%     plot(rsample, K^2 - 2*mu*V(rsample)/hbar^2*unitfactor1- (lsample+1/2)^2./rsample.^2>0);
% end

r = linspace(rmin,rmax, 1000);


for i=1:length(k_array)
    k = k_array(i);
fsum = 0;

r0_array = [];

for ell=0:lmax
    b_r = k^2 - 2*mu*V(r)/hbar^2*unitfactor1- (ell+1/2)^2./r.^2>0;
    
   % if sum(b_r)>0

    %r_local = r(b_r);
   % r0 = r_local(1);
     

    % find outer zeros of the temp expression below
    temp = k^2 - 2*mu*V(r)/hbar^2*unitfactor1- (ell+1/2).^2./r.^2;
    
    % check where the expression is positive
    b_r =temp>0;
    
    % edge detect by shifting by one index
    b_r_shift = circshift(b_r, 1, 2);
    b_r_shift(:,1) = 0;
    edge = ((b_r) & (~b_r_shift)) ;
    
    % detect first zero in flipped edge array ( = outer zero)
    [row, col] = find(flip(edge,2),length(ell), "first");
    % assign r0
    r0 = r( -col + length(r) + 1)';

    


    b_r = r>=r0;
    r0_array = [r0_array,r0];
    r_local = r(b_r);
    
    if ell==6
        diagnostic.temp = temp;
        diagnostic.b_r = b_r;
        
    end

    





    eta = (ell + 1/2)*pi/2 - k*r0 ...
    * trapz(r_local, ...
    sqrt(  k.^2 - 2*mu*V(r_local)/hbar^2 * unitfactor1 - (ell+1/2)^2./r_local.^2) - k );

    fpart = (2*ell+1)*(exp(2*1j*eta) - 1)*legendreP(ell, cos(theta));   
    fsum = fsum + fpart;
   % else
%
       % warning(['domain too small, no contribution from ', sprintf('l = %.0f', ell)]);
   % end

    if r0==rmin
        warning('r0 might be wrong');
        plot(r, k^2 - 2*mu*V(r)/hbar^2*unitfactor1- (ell+1/2)^2./r.^2)
        plot(r, V(r));
    end
    disp(r0);
end
res(i, :) = fsum/(2*1j*k);

disp(k);

end

diagnostic.r0 = r0_array;




end

