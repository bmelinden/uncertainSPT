function W=hmmInit_P_dat(tau,R,Ddt_init,A_init,p0_init,dat)
% W=EMhmmInit_P_dat(tau,R,Ddt_init,A_init,p0_init,dat)
%
% Initialize a diffusive HMM model with 
% tau,R    : blur parameters
% Ddt_init : diffusion constant*timestep
% A_init   : transition matrix 
% p0_init  : initial state probability
% dat      : trajectory data, from EMhmm.preprocess
%
% Number of hidden states given by the length of Ddt_init.
% ML 2016-07-04

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMhmm.hmmInit_P_dat, initialize variational diffusive HMM
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


W=struct;   % model struct
W.N=length(Ddt_init);
W.dim=dat.dim;

W.one  = dat.one;
W.end  = dat.end+1;
W.lnL=0;    % log likelihood
W.tau=tau;
W.R=R;
W.pOcc=zeros(1,W.N);

% parameter subfield
W.P=struct; 
W.P.lambda=2*Ddt_init;
W.P.A=A_init;
W.P.p0=p0_init;

% hidden path subfield, with no Infs or NaNs
W.Y=struct;
W.Y.mu   = dat.x;   % <y(t)>
W.Y.mu(W.end,:)=dat.x(W.end-1,:); % unobsered last position

% fill out missing positions by linear interpolation
ind0=find( isfinite(dat.v(:,1)));
ind1=find(~isfinite(dat.v(:,1)));
for d=1:W.dim
    W.Y.mu(ind1,d)=interp1(ind0,W.Y.mu(ind0,d),ind1,'linear','extrap');
end

W.Y.sig0 = dat.v;   % var(y(t))
W.Y.sig0(W.end,:) = dat.v(W.end-1,:);   % var(y(t))
W.Y.sig0(dat.v==inf)=mean(dat.v(isfinite(dat.v).*(dat.v>0)>0));
if(~isempty(find(~isfinite(W.Y.sig0(:)),1)))
    warning('initial guess generated non-finite position variances, possibly due to dat.v=inf at the end or beginning of trajectories.')
    W.Y.sig0(W.Y.sig0==inf)=mean(dat.v(isfinite(dat.v).*(dat.v>0)>0));
    disp('removed non-finite W.Y.sig0 elements')
end
W.Y.sig1 = zeros(size(dat.v)); % cov(y(t),y(t+1))        
W.Y.MeanLnQ=0;     % entropy of q(y) = <ln q(y)>

% initialize hidden state field
W.S=struct;
W.S.pst=ones(size(dat.x,1),W.N)/W.N;
W.S.wA=ones(W.N,W.N);



