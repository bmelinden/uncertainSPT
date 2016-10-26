function [x,cdf]=empiric_cdf(x)
% [x,cdf]=empiric_cdf(x)
% compute the empirical cumulative distribution function of the
% observations x.

x=reshape(sort(x(:)),1,numel(x));
cdf=(1:numel(x))/numel(x);

