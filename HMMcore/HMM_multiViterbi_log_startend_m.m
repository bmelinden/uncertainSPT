% function S=HMM_multiViterbi_log_m(lnQ,lnH,iStart,iEnd) 
% most likely trajectories by the Viterbi algorithm, using a one-array
% representation of many trajectories, and the log of transition
% matrix lnQ and emission likelihood lnH.  iStart and iEnd are index lists
% of start- and end-points of all single trajectories (using Matlabs index
% conventions, i.e., with first elements indexed by 1). This is the matlab
% version of HMM_multiViterbi_log_startend.c, created for debugging
% purposes, but very much slower.
%
% ML 2015-12-15

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rHMM_multiViterbi_log_startend_m.m, part of HMMcore/
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
%% start of actual code
function S=HMM_multiViterbi_log_startend_m(lnQ,lnH,iStart,iEnd) 
[T,N]=size(lnH);
% algorithm version without large lnpt array
lnPP=zeros(1,N);
lnP0=zeros(1,N);
lnP1=zeros(1,N);
MaxPrev=zeros(T,N,'uint16'); % state variable
S=zeros(T,1,'uint16');

for n=1:length(iEnd)
    tStart=iStart(n);
    tEnd=iEnd(n);
    
    lnP1=lnH(tStart,:)-mean(lnH(tStart,:)); % initial probability, not normalized
    for tV=(tStart+1):tEnd
        lnP0=lnP1;
        lnP1=zeros(1,N);        
        for jV=1:N      % current target state
            for kV=1:N  % previous state
                % probability of most likely path that ends with kV -> jV
                lnPP(kV)=lnP0(kV)+lnQ(kV,jV)+lnH(tV,jV); 
            end
            % probability of previous state before ending up at jV.
            [lnP1(jV),         MaxPrev(tV,jV)]=max(lnPP);
        end
        lnP1=lnP1-mean(lnP1); % rescale to avoid numerical under- or overflow.
    end

    [~,S(tEnd)]=max(lnP1);    
    for tV=tEnd-1:-1:tStart
        S(tV)=MaxPrev(tV+1,S(tV+1));
    end
end


