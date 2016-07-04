% [lnZ,wA,pst]=HMM_multiForwardBackward_startend_m(Q,H,iStart,iEnd,doBackward) 
%
% performs forward-backward sweeps on HMM-type time series stacked on top
% of each other, as described by start- and end indices.
%
% Q    : transition matrix (not necessarily normalized) transition matrix
% H    : emission likelihood, including initial state probabilities on
%        appropriate rows 
% iStart,iEnd : indices to individual time series, e.g., individual
%        trajectories run from  iStart(1) -> iEnd(1), iStart(2) -> iEnd(2),
%        etc. This is so that the algorithm know at what lines the
%        transition matrix should be omitted from the calculations.
%        (NOTE: entries not included in the index ranges are ignored.)

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
% M.L. 2015-12-15

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rHMM_multiForwardBackward_m.m, part of HMMcore/
% =========================================================================
% 
% Copyright (C) 2013 Martin Lindén, E-mail: bmelinden@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% Additional permission under GNU GPL version 3 section 7
% 
% If you modify this Program, or any covered work, by linking or combining it
% with Matlab or any Matlab toolbox, the licensors of this Program grant you 
% additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
