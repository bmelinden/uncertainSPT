function [x,y,Q,yTick,qTick]=QQplot_N01(x)
% [x,y,Q,yTick,qTick]=QQplot_N01(x)
% A Q-Q probability plot of the data x against a N(0,1) distribution.
% 

% empirical quantiles
x=sort(x);                  % empirical quantiles of x
Q=(1:numel(x))/numel(x);    % empirical CDF of x

y=norminv(Q,0,1);           % corresponding quantiles from a N(0,1)

qTick=[1e-3 0.01 0.05 0.1 0.25 0.5 0.75 0.9 0.95 0.99 0.999];
qTick=[1e-3 0.05 0.5 0.95 0.999];
qTick=qTick(qTick>min(Q) & qTick<max(Q));
yTick=norminv(qTick,0,1);

