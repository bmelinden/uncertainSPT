function [BS,ind,BS2] = parameterBootstrap(W,dat,Nbs,dt,display,varargin)
% [BS,ind,BS2] = EMhmm.parameterBootstrap(W,dat,Nbs,dt,display,args)
% Parameter bootstrap 
%
% W     : input HMM model (recommend converged model)
% dat   : data to resample
% Nbs   : number of bootstrap resamples
% dt    : time step for parameter estimate (optional: default 1);
% display : if true, announce convergence of each replica (default: false)
% args  : additional arguments are passed to the parameterEstimate
%         function, e.g., '2state',Dthr, to compute parameters of an
%         aggregated -state model.
%
% BS    : bootstrap struct with bootstrap parameter estiamtes. Last index
%         is the bootstrap iteration. See EMhmm.parameterEstimate for the
%         meaning of the fields.
% ind   : bootstrap indices, one permutation per row.
% BS2   : another bootstrap struct, with parameters from the 2-state
%         collapsed model using the diffusion threshold (only if
%         '2state',Dthr is given).
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

if(~exist('dt','var'))
    dt=1;
end
if(~exist('display','var'))
    display=false;
end

Pargs={};
if(nargin>5) % the additional arguments are present
    Pargs=varargin;
end
[BS,BS2]=EMhmm.parameterEstimate(W,dt,Pargs{:});
% initialize
f=fieldnames(BS);
ind=zeros(Nbs,length(W.i0));
for v=1:length(f)
    [rows,cols]=size(BS.(f{v}));
    BS.(f{v})=zeros(rows,cols,Nbs);
end
P=cell(1,Nbs);
if(isempty(BS2))
    do2state=false;
else
    do2state=true;
    f2=fieldnames(BS2);
    for v=1:length(f2)
        [rows,cols]=size(BS2.(f2{v}));
        BS2.(f2{v})=zeros(rows,cols,Nbs);
    end
    P2=cell(1,Nbs);
end

disp(['Converging ' int2str(Nbs) ' bootstrap replicas.'])
parfor k=1:Nbs
%for k=1:Nbs
    t0=tic;
    % resample and recoverge
    [Xb,Wb,trj]=EMhmm.reorder_trj(dat,W);
    Wb=EMhmm.MLEconverge(Wb,Xb,'display',0);
    if(do2state)
        [P{k},P2{k}]=EMhmm.parameterEstimate(Wb,dt,Pargs{:});
    else
        P{k}=EMhmm.parameterEstimate(Wb,dt,Pargs{:});
    end
    ind(k,:)=trj;
    if(display)
        disp(['Converged replica ' int2str(k) ' in ' num2str(toc(t0)) ' s.'])
    end
end
disp('Done.')
for k=1:Nbs
    for v=1:length(f)
        BS.(f{v})(:,:,k)=P{k}.(f{v});
    end
    if(do2state)
        for v=1:length(f2)
            BS2.(f2{v})(:,:,k)=P2{k}.(f2{v});
        end
    end
end
