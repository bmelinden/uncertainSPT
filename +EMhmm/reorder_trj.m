function [Xb,Wb,trj]=reorder_trj(X0,W0,trj)
% [Xb,Wb,trj]=reorder_trj(X0,W0,trj)
% Assemble a bootstrapped data-model pair from a converged RMhmm model,
% conserving the current model parameters and variational distributions.
%
%
% X0    : original EMhmm data struct
% W0    : original EMhmm model, converged with X0
% trj   : (optional) trajectory reordering indices. If empty or not given,
%          a random resampling with replacement of the same number of
%          trajectories as X0 is generated, suitable for bootstraping. 
% 
% Xb    : re-ordered trajectory struct
% Wb    : re-ordered EMhmm model
% trj   : trajectory indices for reordering
%
% ML 2016-08-19

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reorder_trj, trajectory reordering for diffusive HMM analysis
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
Ntrj=length(W0.i0); % number of trajectories in the model

% resampled trajectory indices if not given
if(~exist('trj','var') || isempty(trj))
    trj=sort(ceil(Ntrj*rand(1,Ntrj)));  
end
T=X0.T(trj);
dim=X0.dim;

% rebuild trajectories
doMisc=isfield(X0,'misc');

Xb=struct;
Xb.dim=dim;
Xb.T=T;
Xb.x=zeros(sum(T+1),dim);
Xb.v=zeros(sum(T+1),dim);
Xb.i0=zeros(1,length(T),'double');
Xb.i1 =zeros(1,length(T),'double');

ind=1;
for k=1:length(trj)
    k0=trj(k);
    
    x=X0.x(X0.i0(k0):X0.i1(k0),1:dim);
    v=X0.v(X0.i0(k0):X0.i1(k0),1:dim);
    Tx=size(x,1);

    Xb.i0(k)=ind;
    Xb.i1(k)=ind+Tx-1;
    ind=ind+Tx;
    Xb.x(Xb.i0(k):Xb.i1(k),1:dim)=x; 
    Xb.v(Xb.i0(k):Xb.i1(k),1:dim)=v;
    if(doMisc)
       Xb.misc( Xb.i0(k):Xb.i1(k),:)=X0.misc( X0.i0(k0):X0.i1(k0),:);
    end
    ind=ind+1;
end

% rebuild model in the same way
% first initialize
Wb=EMhmm.init_P_dat(W0.tau,W0.R,W0.P.lambda/2,W0.P.A,W0.P.p0,Xb);
% then transfer variational fields
Wb.S.wA=W0.S.wA; % requires no permutation
for k=1:length(trj)
    k0=trj(k);
    
    Wb.Y.mu(Wb.i0(k):Wb.i1(k),:)  =W0.Y.mu(W0.i0(k0):W0.i1(k0),:);
    Wb.Y.sig0(Wb.i0(k):Wb.i1(k),:)=W0.Y.sig0(W0.i0(k0):W0.i1(k0),:);
    Wb.Y.sig1(Wb.i0(k):Wb.i1(k),:)=W0.Y.sig1(W0.i0(k0):W0.i1(k0),:);
    Wb.S.pst(Wb.i0(k):Wb.i1(k),:) =W0.S.pst(W0.i0(k0):W0.i1(k0),:);
end
Wb.pOcc=rowNormalize(sum(Wb.S.pst,1));
