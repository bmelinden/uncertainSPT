% function S=HMM_multiViterbi_log_m(lnQ,lnH,iStart,iEnd) 
% most likely trajectories by the Viterbi algorithm, using a one-array
% representation of many trajectories, and the log of transition
% matrix lnQ and emission likelihood lnH. iStart and iEnd are index lists
% of start- and end-points of all single trajectories (using Matlabs index
% conventions, i.e., with first elements indexed by 1).
%
% ML 2015-12-15
