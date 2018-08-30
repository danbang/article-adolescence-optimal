function slope = quickSlope(delta,choice)
% function slope = quickSlope(delta,choice)
% computes sensitivity as specified in Bahrami et al (2010) Optimally interacting minds 
dval         = delta;
dsteps       = unique(dval);
for k = 1:numel(dsteps)
    j = dsteps(k);
    fsM(k) = length(choice(dval==j & choice==1));
    fsN(k) = length(choice(dval==j));
end
y                = fsM'./fsN';
bhat             = glmfit(dsteps,[y ones(size(y))],'binomial','link','probit');
slope_sd         = 1/bhat(2);
slope            = 1 ./ (sqrt(2*pi) .* slope_sd);
end