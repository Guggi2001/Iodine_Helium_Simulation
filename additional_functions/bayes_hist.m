function [h, sigma_h, centers, binwidth, barplot, errorplot] = bayes_hist(plot_samples, xinterval, no_plot, col, edges)

if nargin<2
    xinterval = [min(plot_samples), max(plot_samples)];
end

if nargin<3
    no_plot = false;
    col = 'b';
end

xmin = xinterval(1);
xmax = xinterval(2);
N = numel(plot_samples);

if nargin<5
[n, edges] = histcounts(plot_samples, 'BinLimits', [xmin xmax]);
else
[n, edges] = histcounts(plot_samples, edges);
end
nbins = length(n);
binwidth = edges(2)-edges(1);
centers = edges(1:end-1)+ binwidth/2;




p = (n+1)/(N+nbins + 1);
sigma_p = sqrt( (n+2)/(N + nbins + 2).*p - p.^2 );


h = p/binwidth;
sigma_h = sigma_p/binwidth;


if ~no_plot



errorplot = errorbar(centers, h,sigma_h,'CapSize', 2,'Color',col,'LineWidth',1.2,'Marker','none', 'LineStyle','none', 'HandleVisibility', 'off') ; 
 
hold on;
barplot = bar(centers,h, 1, 'FaceAlpha',0.0,'FaceColor', col,'EdgeColor', col, 'EdgeAlpha',0.8);
else
    barplot = nan;
    errorplot = nan;
end


end

