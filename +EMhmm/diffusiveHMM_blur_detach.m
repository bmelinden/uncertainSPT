function [x,v,s,y,ym]=VB6_diffusiveHMM_blur_detach(p0,A,pE,Ddt,tE,locErrG,dim,T,numTrj,pM)
% [x,v,s,y]=VB6_diffusiveHMM_blur_detach(p0,A,pE,Ddt,tE,locErrG,dim,T,numTrj,pM)
%
% Simulate diffusive HMM data. 
% x(t,:)  : measured positions
% v(t,:)  : variance of localization errors in x
% s(t)    : hidden state sequence
% y(t,:)  : true positions
% ym(t,:) : motion-averaged true positions
% Motional blur and localization errors are added according to the
% Michalet&Berglund model with independent Gauyssian noise, whose std are
% iid random variables.
% Cell vectors are returned if more than one trajectory is simulated.
%
% p0 : initial state distribution. Default: stationary state of A. p0 is
%      automatically normalized.
% A  : transition matrix, using the convention 
%      A(i,j)=p(s_t=j|s_{t-1}=i), so the evolution of the hidden states is
%      given by p(t)=p(t-1)*A, where p(t) is a row vector
% pE : trajectory end probability vector, pE(j)=p(t is last point|s_t=j).
%      Default 0, and a scalar is interpreted as same for all states.
%      Minimum trajectory length is 2.
% Ddt: diffusion*timestep for the various hidden states, Ddt(s(t)==j)=D(j)*dt
% tE : exposure time, in units of the sample time dt. This code simulates
%      uniform illumination, and so R=tE/6, tau=tE/2 is implied. tE=0 means
%      no blur, only localization errors.
% locErrG=[Em Estd] : parameterization of static localization error
%      distribution. Here, each localization error in each component are
%      independent Gaussians with zero mean and iid standard devitions St
%      which are gamma-distributed with mean Em and standard deviation
%      Estd, both >0.
% dim: spatial dimension of the output data (default 2)
% T  : maximum trajectory length(s). Only the first numTrj entries are
%      used, and if numTrj>length(T), then the entries are cycled. Default: 100. 
% numTrj : number of trajectories to simulate. Default: 1.
% pM : fraction of missing positions (x(t)=NaN, v(t)=Inf). Default : 0.
%
% ML 2015-11-30, bmelinden@gmail.com
% ML 2015-12-09 :  missing data points

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMhmm.diffusiveHMM_blur_detach.m, simulate d-dimensional diffusion with 
% multiple diffusion constants and random errors
% =========================================================================
% 
% Copyright (C) 2015 Martin Lindén
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

%% parameter handling
[N,NN]=size(A);
if(N~=NN)
    error('VB6_diffusiveHMM_blur_detach requires square transition matrix')
end
clear NN

if( isempty(p0) || length(p0)==1 )
    A0=A^10000;
    p0=A0(1,:);
end
p0=p0/sum(p0);

if(isempty(pE))
    pE=zeros(1,N);
elseif(length(pE)==1)
    pE=pE*ones(1,N);
end
if( tE <0 || tE >1)
    error('Exposure time tE must be in the range [0,1].')
end
tau=tE/2;
R=tE/6;

if(length(locErrG)~=2)
    error('locErrG should be specified as [Em Estd].')
end
ErrShape=(locErrG(1)/locErrG(2))^2; % gamma shape parameter
ErrScale=locErrG(2)^2/locErrG(1);   % gamma scale parameter

if(~exist('dim','var') || isempty(dim)); dim=2; end
if(~exist('T','var') || isempty(T)); T=100; end
if(~exist('numTrj','var') || isempty(numTrj)); numTrj=1; end
if(~exist('pM','var') || isempty(pM)); pM=0; end

Ddt=reshape(Ddt,length(Ddt),1);

%% initialization
x=cell(1,numTrj);
v=cell(1,numTrj);
s=cell(1,numTrj);
y=cell(1,numTrj);
%% start simulation
NT=length(T);
cumA=cumsum(A,2);
cum0=cumsum(p0);
for m=1:numTrj
    % initialize trajectory
    Ttrj=T(1+mod(m-1,NT)); % max length of this trajectory
    if(Ttrj<inf)
        S=zeros(Ttrj,1);
    else
        S=zeros(10,1);
    end
    % hidden state trajectory
    S(1)=find(rand<cum0,1); % initial state
    for t=2:Ttrj
        ras=rand;
        S(t)=find(ras<cumA(S(t-1),:),1);
        if(rand<pE(S(t))) % then terminate trajectory here
            S=S(1:t);
            break
        end
    end
    Ttrj=length(S);
    
    % pure diffusion-trace
    dY=zeros(Ttrj,dim);
    for k=1:dim
        dY(:,k)=randn(Ttrj,1).*sqrt(2*Ddt(S));
    end
    Y=[zeros(1,dim); cumsum(dY,1)];
    
    % iid localization errors
    dX=gamrnd(ErrShape,ErrScale,size(Y,1)-1,dim);
    
    % add motional blur     
    YM=Y(1:Ttrj,:)*(1-tau)+Y(2:end,:)*tau;
    for k=1:dim
        YM(:,k)=YM(:,k)+randn(Ttrj,1).*...
            sqrt(2*Ddt(S)*(tau*(1-tau)-R));
    end
    % and localization noise
    for k=1:dim
        X(:,k)=YM(:,k)+randn(Ttrj,1).*dX(:,k);
    end
    
    M=rand(size(X))<pM;
    % insert missing positions
    if(pM>0)
        M=rand(size(X,1),1)<=pM;
        X(M,:) =0;
        dX(M,:)=Inf;
    end
    
   % if(numTrj==1)
    %    x=X;s=S;y=Y;v=dX.^2;
    %else
        x{m}=X;s{m}=S;y{m}=Y;ym{m}=YM;v{m}=dX.^2;
    %end
end
