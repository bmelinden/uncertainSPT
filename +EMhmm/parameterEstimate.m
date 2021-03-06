function [est,est2]=parameterEstimate(W,dt,varargin)
% [est,est2]=parameterEstimate(W,dt,...)
% Estimate some model properties and parameters.
%
% minimal operation:
%
% W     : converge HMM model struct
% dt    : time step of the data (default: 1)
%
% est   : parameter struct with fields
%             D: diffusion constant
%        lambda: 2*D*dt
%            p0: initial state probability
%          pOcc: total occupancy (from variational state distribution)
%           pSS: steady state of the transition matrix A
%             A: transition matrix
%    dwellSteps: mean dwell time in units of time step
%     dwellTime: mean dwell times in units of time
% est2  : same as est, but for the 2-state coarse-grained model described
%         below.
%
% optional input: parameterEstimate(W,dt,'2state',Dthr) produces additional
%               parameter estimates est2 for a coarse-grained model with a
%               slow (D<=Dthr) and a fast (D>Dthr) state, by simply adding
%               up the summary statistics into two groups. If either of the
%               two groups are empty, all parameter estimates are NaN.
%
% ML 2017-01-27

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameterEstimate, estimate some model properties from diffusive HMM
% =========================================================================
% 
% Copyright (C) 2017 Martin Lindén
% 
% E-mail: bmelinden@gmail.com
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
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by linking or combining it
%  with Matlab or any Matlab toolbox, the licensors of this Program grant you 
%  additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.

if(~exist('dt','var'))
    dt=1;
end

k=0;
slowFastAggregate=false;
while(k<nargin-2)
    k=k+1;
    if(strcmp(varargin{k},'2state'))
       k=k+1;
       Dthr=varargin{k};
       slowFastAggregate=true;
    else
        error(['Argument ' varargin{k} ' not recognized.'])
    end
end

est=struct;
est.D = W.P.lambda/2/dt;
est.lambda=W.P.lambda;
est.p0= W.P.p0;
est.pOcc=rowNormalize(sum(W.S.pst,1));
est.A = W.P.A;
[est.pSS,allOK] =EMhmm.steadyStateFromA(est.A);
est.dwellSteps= 1./(1-diag(est.A)'); % mean dwell times [steps]
est.dwellTime = dt./(1-diag(est.A)'); % mean dwell times [time units]

Tstates=sum(W.S.pst(:));
est.arrRate=sum(W.S.wA.*(1-eye(W.N)),1)/Tstates;
est.depRate=sum(W.S.wA.*(1-eye(W.N)),2)'/Tstates;

if(~allOK)
   warning('Problem computing full steady state. wA =') 
   disp(num2str(W.S.wA))  
   disp('A=')
   disp(num2str(W.P.A))
end

% print and plot MLE results with thresholds
est2=struct;
if(slowFastAggregate)
    iSlow=find(est.D<=Dthr);
    iFast=find(est.D>Dthr);
    est2.isSlow=est.D<=Dthr;
    if(isempty(iSlow) || isempty(iFast) )
        warning('EMhmm.parameterEstimate 2-state coarsegraining: D threshold outside D interval')        
        est2.D=nan(1,2);
        est2.p0=nan(1,2);
        est2.pOcc=nan(1,2);
        est2.A=nan(2,2);
        est2.pSS=nan(1,2);
        est2.dwellSteps=nan(1,2);
        est2.dwellTime=nan(1,2);
        est2.arrRate=nan(1,2);
        est2.depRate=nan(1,2);
    else
        est2.D=[est.D(iSlow)*rowNormalize(est.pOcc(iSlow))'  est.D(iFast)*rowNormalize(est.pOcc(iFast))'];
        est2.p0=[sum(est.p0(iSlow))  sum(est.p0(iFast))  ];
        est2.pOcc=[sum(est.pOcc(iSlow))  sum(est.pOcc(iFast))  ];
        wA=zeros(2,2);
        wA(1,1)=sum(sum(W.S.wA(iSlow,iSlow)));
        wA(1,2)=sum(sum(W.S.wA(iSlow,iFast)));
        wA(2,1)=sum(sum(W.S.wA(iFast,iSlow)));
        wA(2,2)=sum(sum(W.S.wA(iFast,iFast)));
        est2.A=rowNormalize(wA);
        [est2.pSS,allOK]=EMhmm.steadyStateFromA(est2.A);
        est2.dwellSteps= 1./(1-diag(est2.A)'); % mean dwell times [steps]
        est2.dwellTime = dt./(1-diag(est2.A)'); % mean dwell times [time units]
        est2.arrRate=sum(wA.*(1-eye(2)),1)/Tstates;
        est2.depRate=sum(wA.*(1-eye(2)),2)'/Tstates;

        if(~allOK)
            warning('Problem computing reduced steady state. wA = ')
            disp(num2str(wA))
            disp('A=')
            disp(num2str(W.P.A))
        end
    end
end
