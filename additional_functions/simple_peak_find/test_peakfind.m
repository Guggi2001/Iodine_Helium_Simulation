x = -10:0.1:10;

y = g(4, 2, x) + g(-5,1,x) - g(-2, 6, x);

y = y.*(1 + (rand(size(y))-0.5)/2) + rand(size(y))*1;

y = movmean(y,10,2)

plot(x,y)


ids = findLocalMaxima(y,31,90,60);

hold on
scatter(x(ids), y(ids));
