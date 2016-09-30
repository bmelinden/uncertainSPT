function W=MLEparameterUpdate(W,dat)
% one round of EM iterations in a diffusive HMM, including handling missing
% data, and adding conjgate prior fields
% P0.wA, P0.w0, P0.n, P0.c if the model W has not prior field P0.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMhmm.MLEparameterUpdate, variational diffusive HMM, postion estimator
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

%% parameter update
W.P.A=rowNormalize(W.S.wA);
W.P.p0=rowNormalize(sum(W.S.pst(W.one,:),1));

% check for unoccupied rows
wAemptyRows=find((sum(W.S.wA,2)==0))';
if(~isempty(wAemptyRows)) % a very non-invasive regularization, no new transitions
   W.P.A=rowNormalize(W.S.wA+10*eps*eye(W.N));
   W.P.p0=rowNormalize(sum(W.S.pst(W.one,:),1)+10*eps);
   %warning(['State(s) ' int2str(wAemptyRows) ' unoccupied, MLEparameterUpdate adding 10*eps pseudocounts to avoid NaNs'])
end


% second attempt at diffusion constants
MLEopt = optimoptions('fminunc','GradObj','off','TolX',1e-14,...
    'MaxIter',1e5,'TolFun',1e-14,'Algorithm','quasi-newton',...
    'MaxFunEvals',1e10,'DerivativeCheck','off','Display','none');

indS=1:sum(dat.T+1);
indS(W.end)=0;
indS=indS(indS>0);
dy2=sum((W.Y.mu(indS+1,:)-W.Y.mu(indS,:)).^2 ...
    +W.Y.sig0(indS,:)+W.Y.sig0(indS+1,:)-2*W.Y.sig1(indS,:),2);
dxy2=(dat.x(indS,:)-W.Y.mu(indS,:)*(1-tau)-W.Y.mu(indS+1,:)*tau).^2 ...
        +(1-tau)^2*W.Y.sig0(indS,:)...
        +tau^2*W.Y.sig0(indS+1,:)...
        +2*tau*(1-tau)*W.Y.sig1(indS,:);
dxy2(~isfinite(dat.v(indS,:)))=0; %%% ad-hoc warning, use pO instead
for j=1:W.N
    pt=W.S.pst(indS,j);
    nj=W.dim/2*sum(pt);    
    cj=0.5*sum(pt.*dy2);
    
    if(nj>0) % only MLE update if there is any occupancy on which to base the update
        Fj=@(L)(-nj*log(L)-cj/L-0.5*sum(pt.*sum(log1p(beta*L./dat.v(indS,:))+dxy2./(dat.v(indS,:)+beta*L),2)));
    
        mLogLj=@(logLam)(-Fj(exp(logLam)));
        W.P.lambda(j)=exp(fminunc(mLogLj,log(W.P.lambda(j)),MLEopt));
    else
        W.P.lambda(j)=1e100; % set rediculously large value to put unoccupied state last in ordered models
    end
end



