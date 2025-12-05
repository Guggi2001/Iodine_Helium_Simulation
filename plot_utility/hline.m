function hl  = hline(y,color, z, linestyle)

%HLINE
if nargin<3
    z = 0;
end

if nargin<2
    color = [0.5, 0.5, 0.5];
end

if nargin<4
    linestyle = '--';
end



xlimits = get(gca, 'Xlim');

if z==0
    for i = y
        hl = plot(xlimits, [i,i],linestyle, 'color', color, 'linewidth',1.5, 'HandleVisibility', 'off');
    end
else
    for i = y
        hl = plot3(xlimits,[i,i],[z,z],linestyle, 'color', color,'linewidth', 1.5, 'HandleVisibility', 'off');
    end
end
end

