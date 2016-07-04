function [W,sMaxP,sVit,WS]=hiddenStateUpdate(W,dat)
% [W,sMaxP,sVit,WS]=hiddenStateUpdate(W,dat)
% one round of hidden state EM iteration in a diffusive HMM, with possibly
% missing position data
%
% sMAxP and sVit are only computed if asked for by output arguments.
% WS : struct containing the whole workspace at the end of the function
% call. Expensive, computed only when asked for.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMhmm.hiddenStateUpdate, variational state update for diffusive HMM
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
%% assemble point-wise weights
lnH=zeros(size(W.S.pst ));   
for t=1:length(W.one)
    % initial state distribution
    lnH(W.one(t),:)=log(W.P.p0);
    % -d/2*log(lambda_j)
    indY=W.one(t):W.end(t); % indices to all hidden positions
    indS=indY(1:end-1);     % indices to hidden states and observed positions
    indO=indS(isfinite(dat.v(indS,1))); % indices to non-missing observed data
    TO=numel(indO); % number of observed positions in this trj
    for j=1:W.N
        lnH(indS,j)=lnH(indS,j)-W.dim/2*log(W.P.lambda(j));
    end
    % average steplength variances
    sumK_dy2= sum(diff(W.Y.mu(indY,:)).^2 ...
        +W.Y.sig0(indS,:) ...
        +W.Y.sig0(1+indS,:) ...
        -2*W.Y.sig1(indS,:),2);
    lnH(indS,:)=lnH(indS,:)-sumK_dy2/2*(1./W.P.lambda);
    % dynamic localization error
    %%%%% got this far in missing point correction
    %dx2tk=(dat.x(indS,:)-W.Y.mu(indS,:)*(1-tau)-W.Y.mu(1+indS,:)*tau).^2 ...
    %    +(1-tau)^2*W.Y.sig0(indS,:)+tau^2*W.Y.sig0(1+indS,:)...
    %    +2*tau*(1-tau)*W.Y.sig1(indS,:);
    % version with potentially missing points
    %   for j=1:W.N
    %    vtk_blj=dat.v(indS,:)+beta*ones(dat.T(t),1)*W.P.lambda(j)*ones(1,W.dim);
    %    lnH(indS,j)=lnH(indS,j)-0.5*sum(log(vtk_blj)+dx2tk./vtk_blj,2);
    %end
    dx2tk=(dat.x(indO,:)-W.Y.mu(indO,:)*(1-tau)-W.Y.mu(1+indO,:)*tau).^2 ...
        +(1-tau)^2*W.Y.sig0(indO,:)+tau^2*W.Y.sig0(1+indO,:)...
        +2*tau*(1-tau)*W.Y.sig1(indO,:);
    % dynamic localization variance
    for j=1:W.N
        vtk_blj=dat.v(indO,:)+beta*W.P.lambda(j)*ones(TO,W.dim);
        %lnH(indO,j)=lnH(indO,j)-0.5*sum(log(vtk_blj)+dx2tk./vtk_blj,2);
        blvm1j=beta*W.P.lambda(j)./dat.v(indO,:); 
        lnH(indO,j)=lnH(indO,j)-0.5*sum(log1p(blvm1j)+dx2tk./vtk_blj,2);
        if(~isempty(find(~isfinite(lnH(:,j)),1)))
           error('NaN is lnH!') 
        end
    end
    % using the log(1+...) notation, might be numerically better since it
    % avoids adding constants that are then just subtracted off again in
    % lnHmax. However, log1p might be slower than just log...
end
lnHmax=max(lnH,[],2);
lnH=lnH-lnHmax*ones(1,W.N);
H=exp(lnH);
H(W.end,:)=0;
%% forward-backward iteration
%[ln3,wA3,ps3]=HMM_multiForwardBackward_g1(W.P.A,H,dat.end);
[lnZ,W.S.wA,W.S.pst]=HMM_multiForwardBackward_startend(W.P.A,H,dat.one,dat.end);
W.S.lnZ=lnZ+sum(lnHmax);
W.pOcc=rowNormalize(sum(W.S.pst,1));
%% likelihood lower bound after s
W.lnL=W.S.lnZ-W.Y.MeanLnQ;
%% path estimates
if(nargout>=2) % compute sequence of most likely states
   [~,sMaxP]=max(W.S.pst,[],2); 
   sMaxP(W.end)=0;
end
if(nargout>=3) % compute Viterbi path, with a small offset to avoid infinities
   sVit=HMM_multiViterbi_log_startend(log(W.P.A+1e-50),log(H+1e-50),W.one,W.end-1);
end
if(nargout>=4)
   fname=['foo_' int2str(ceil(1e5*rand)) '.mat'];
   save(fname);
   WS=load(fname);
   delete(fname);
end
