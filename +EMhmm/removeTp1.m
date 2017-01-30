function [t,X]=removeTp1(dat,X)
% insert NaNs in the T+1-gaps from a vector representation X, and
% compute corresponding time step indices. A simple way to plot
% trajectory data in a single command and still get indications of
% where trajectories start and end.
%
%
t=(1:dat.i1(end)+1)';
X(dat.i1+1,:)=NaN;

