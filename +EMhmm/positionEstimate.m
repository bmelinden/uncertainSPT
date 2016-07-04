function xEst=positionEstimate(W,dat)
% xEst=positionEstimate(W,dat)
% Estimate actual positions (posterior mean position) based on input data
% and converge HMM model.
%
% ML 2016-06-07

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMhmm.positionEstimate, variational diffusive HMM, postion estimator
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
%% start of actual code


beta=W.tau*(1-W.tau)-W.R;
tau=W.tau;

bvl=zeros(size(W.S.pst,1),1); % beta*lambda_t/v_t
ym=zeros(size(W.S.pst,1),1);  % y_t*(1-tau)+y_{t+1}*tau
xEst=zeros(size(W.S.pst,1),W.dim); 

% v 2
for d=1:W.dim
   ym(1:end-1,1)=(1-tau)*W.Y.mu(1:end-1,d)+tau*W.Y.mu(2:end,d);
   bvl=(beta./dat.v(:,d))*W.P.lambda; 
   xEst(:,d)=ym.*sum(W.S.pst./(1+bvl),2)...
            +dat.x(:,d).*sum(W.S.pst.*bvl./(1+bvl),2);
end
xEst(W.end,:)=0;               % make sure there is no position estimate for t=T+1

% variance of position estimate
if(0)
vxEst=zeros(size(xEst));
for d=1:W.dim
   
    % conditional variance
    vzt=(dat.v(:,d)*beta*W.P.lambda)./(...
        ones(size(dat.v(:,d)))*beta*W.P.lambda...
        +dat.v(:,1)*ones(size(W.P.lambda)));
    
    mvzt=mean(vzt.*W.S.pst,2);

    % contribution from hidden position uncertainty and states
    varYm=W.Y.sig0(1:end-1,d)*(1-tau)^2+W.Y.sig0(2:end,d)*tau^2 ...
    +2*tau*(1-tau)*W.Y.sig1(1:end-1,d);


    % contribution from hidden states
    

end

end
