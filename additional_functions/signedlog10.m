% signed log functin from forum post
function out = signedlog10(in)
    % naive signed-log
    %out = sign(in).*log10(abs(in));
    
    % modified continuous signed-log
    % see Measurement Science and Technology (Webber, 2012)
    C = 1; % controls smallest order of magnitude near zero
    out = sign(in).*(log10(1+abs(in)/(10^C)));
end

