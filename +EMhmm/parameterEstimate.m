function est=parameterEstimate(W,dt)
% est=parameterEstimate(W)
% Estimate some model properties
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
% 
% ML 2016-08-19

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameterEstimate, estimate some model properties from diffusive HMM
% =========================================================================
% 
% Copyright (C) 2016 Martin Lindén
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

% compute steady state
b=eig(W.P.A);
if(sum(abs(b-1)<10*eps)>1)
    warning('Steady state possibly not unique.')
end
bMax=max(b);
b2nd=max(b(b<bMax));
if(abs(bMax-1)<10*eps && ~isempty(b2nd))
    Nss=log(eps)/log(b2nd);
    ASS=(W.P.A/bMax)^Nss; % this deals with the possibility that the largest eigenvalue is not exactly 1.
    pSS=ASS(1,:);
else
    pSS=NaN(1,W.N);
    warning('Steady state not found.')
end

est.D = W.P.lambda/2/dt;
est.lambda=W.P.lambda;
est.p0= W.P.p0;
est.pOcc=rowNormalize(sum(W.S.pst,1));
est.pSS =pSS;
est.A = W.P.A;
est.dwellSteps= 1./(1-diag(est.A)'); % mean dwell times [steps]
est.dwellTime = dt./(1-diag(est.A)'); % mean dwell times [time units]


