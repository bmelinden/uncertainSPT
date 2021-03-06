function [W,WS]=diffusionPathUpdate(W,dat)
% W=diffusionPathUpdate(W,dat)
% one round of diffusion path update in a diffusive HMM, with possibly
% missing position data.
%
% WS : optional output: workspace at the end of the function

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMhmm.diffusionPathUpdate.m, variational diffusion path update in 
% diffusive HMM
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

%% hidden path: optimized and validated version
nu=zeros(size(dat.x));
lambda_inv=sum((ones(size(dat.x,1),1)*(1./W.P.lambda)).*W.S.pst,2); % sum_j <s(t)==j>/lambda(j)
alpha=zeros(size(dat.x));
for k=1:W.dim % alpha(t,k)=sum_j <s(t)==j>/(v(t,k)+beta*lambda(j))
    alpha(:,k)=sum(W.S.pst./(dat.v(:,k)*ones(1,W.N)...
        +beta*(ones(size(dat.x,1),1)*W.P.lambda)),2);
    % this works with missing points, since v=inf -> alpha=0, and all
    % x-dependent terms are proportional to alpha
end
alphaX=alpha.*dat.x;
alphaX(~isfinite(dat.v))=0;
%%%% to do: replace these ad-hoc corrections with a global indicator
%%%% variable pO, where pO(t)=1 means existing observation, and
%%%% pO(t)=missing data. This is a precursor to having pO being a
%%%% continuous variable to indicate strange points in general

zeros(size(dat.x));
L0=zeros(size(dat.x)); % Lambda matrix diagonal     L0(t)=Lambda(t,t);
L1=zeros(size(dat.x)); % Lambda matrix off-diagonal L1(t)=Lambda(t,t+1)
for t=1:length(W.i0)
    indY=W.i0(t):W.i1(t); % indices to all hidden positions
    indS=indY(1:end-1);     % indices to hidden states and observed positions

    % dynamic localization variance
    nu(indY(1),:)      =alphaX(indS(1),:)*(1-tau);
    nu(indY(2:end-1),:)=alphaX(indS(2:end),:)*(1-tau)...
                       +alphaX(indS(1:end-1),:)*tau;
    nu(indY(end),:)    =alphaX(indS(end),:)*tau;
    
    for k=1:W.dim
        L0(indY(1),k)      =lambda_inv(indS(1  ))  +(1-tau)^2*alpha(indS(1)    ,k);
        L0(indY(end),k)    =lambda_inv(indS(end))  +    tau^2*alpha(indS(end)  ,k);
        L0(indY(2:end-1),k)=lambda_inv(indS(2:end))+(1-tau)^2*alpha(indS(2:end)  ,k)...
                           +lambda_inv(indS(1:end-1))+  tau^2*alpha(indS(1:end-1),k);
                       
        L1(indS,k) = alpha(indS,k)*tau*(1-tau)-lambda_inv(indS);
        
        T=dat.T(t);
        Lt_elements=[L0(indY,k);L1(indS,k);L1(indS,k)];
        Lt_rows=[1:T+1 2:T+1 1:T  ];
        L_cols=[1:T+1 1:T   2:T+1];
        Lt=sparse(Lt_rows,L_cols,Lt_elements,T+1,T+1,3*(T+1));
        W.Y.mu(indY,k)=Lt\nu(indY,k);
        if( ~isempty(find(~isfinite(W.Y.mu(:,k)),1)))
           error('NaN in the path!')
        end
    end    
end
% covariances
logDet=zeros(1,W.dim);
for k=1:W.dim
    if(1)
        [W.Y.sig0(:,k),W.Y.sig1(:,k),logDet(k)]=triSym_triInv_rescale_trjWise(...
            L0(:,k),L1(:,k),W.i0,W.i1,numel(W.i0));
        % a possible fix when something goes really wrong with the recursive
        % tridiagonal inversion, which is not so good numerically.
        if( ~isempty(find(~isfinite(W.Y.mu(:,k)),1)) || ...
                ~isempty(find(~isfinite(W.Y.sig0(:,k)),1)) || ...
                ~isempty(find(~isfinite(W.Y.sig1(:,k)),1)) || ...
                ~isfinite(logDet(k)) )
            % find the offending trajectories and try to correct them
            llD=zeros(1,numel(W.i0));
            for m=1:numel(W.i0)
                ind=W.i0(m):W.i1(m);
                [ss0,ss1,llD(m)]=triSym_triInv_rescale_trjWise(L0(ind,k),L1(ind,k),1,1+W.i1(m)-W.i0(m),1);
                if(~isfinite(llD(m)))
                    LD=spdiags([ L1(ind,k) L0(ind,k) [0;L1(ind(1:end-1),k)]],-1:1,length(ind),length(ind));
                    llD(m)=sum(log(eig(LD)));
                    
                    LDi=inv(LD);
                    ss0=diag(LDi);
                    ss1=[diag(LDi,1);0];
                    W.Y.sig0(ind,k)=ss0;
                    W.Y.sig1(ind,k)=ss1;
                    
                    warning(['nonfinite hidden path model in trj' int2str(m)])
                end
            end
            logDet(k)=sum(llD);
        end
    else 
        % alternative path update that uses only matlabs built-in routines,
        % possibly giving a more stable (if slower) solution.
        [W.Y.sig0(:,k),W.Y.sig1(:,k),logDet(k)]=triSym_triInv_backslash(...
            L0(:,k),L1(:,k),W.i0,W.i1);
    end
    
    %[sig0(:,k),sig1(:,k),lD(k)]=triSym_triInv_rescale_trjWise_m(L0(:,k),L1(:,k),W.i0,W.i1,numel(W.i0));
end
W.Y.MeanLnQ=-sum(W.dim/2*(dat.T+1)*(1+log(2*pi)))+0.5*sum(logDet); % <lnQ(y)>
if(nargout>=2)
   fname=['foo_' int2str(ceil(1e5*rand)) '.mat'];
   save(fname);
   WS=load(fname);
   delete(fname);
end
