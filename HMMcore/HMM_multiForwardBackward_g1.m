% [lnZ,wA,pst]=HMM_multiForwardBackward_m(Q,H,iEnd,doBackward) 
%
% performs forward-backward sweeps on HMM-type time series stacked on top
% of each other. 
%
% Q    : transition matrix (not necessarily normalized) transition matrix
% H    : emission likelihood, including initial state probabilities on
%        appropriate rows 
% iEnd : indices to ends of individual time series, e.g., individual
%        trajectories run from  1 -> iEnd(1), iEnd(1)+2 -> iEnd(2), etc.
%        This is so that the algorithm know at what lines the transition
%        matrix should be omitted from the calculations.
%        (NOTE: entries at rows iEnd+1 are ignored: this is to accomodate
%        blurry HMMs, where ML thought this convention would lead to simple
%        matrix algebra.) 

% doBackward: flag to activate a backward sweep and compute pst (default
%       false).
%        
% lnZ  : log of normalization constant for the forward sweep
% wA   : transition count matrix (requires backward sweep)
% pst  : pst(t,j)= P(s(t)==j), averaged over all state sequences (some
%        extra computational cost in addition to backward sweep).
%        
% This is the matlab implementation of HMM_multiForwardBackward.c, created
% for debugging purposes (much slower). 
%
% M.L. 2015-12-02
