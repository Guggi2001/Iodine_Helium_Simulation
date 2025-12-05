rn = rand(200000,1);

R = 4;

pdf = @(b)2*b/R^2;

ICDF = @(x) sqrt(x*R^2);


samples = ICDF(rn);

figure
b=0:0.1:R;
plot(b, pdf(b));

hold on
bayes_hist(samples, [0,R], false, 'red');

