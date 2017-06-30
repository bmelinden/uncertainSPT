function [W,sMaxP,sVit]=MLEconverge(W,dat,varargin)
% [W,sMaxP,sVit]=MLEconverge(W,dat)
% Run MLE EM iterations on the diffusive HMM W and data dat, until
% convergence.
%
% Output:
% W     : converged HMM model
% sMaxP : trajectory of most likely states 
% sVit  : most likely trajectory of hidden states (from Viterbi algorithm)
% sMaxP and sVit are a little expensive, and only computed if asked for by
% output arguments.
%
% Input :
% W     : EM6 HMM model struct, e.g., created by EMhmm.init_P_dat
% dat   : EM6 data struct, e.g., from EMhmm.preprocess
% optional arguments in the form 'name', value
% Nwarmup   : number of initial iterations where model parameters are kept
%             constant in order to 'burn in' the states and hidden path.
%             Default 5.
% maxIter   : maximum number of iterations (past Nwarmup). Default 5000.
% lnLrelTol : relative convergence criteria for (lnL(n)-lnL(n-1))/|lnL(n)|.
%             Default 1e-8;
% parTol    : convergence criteria for parameters lambda (relative), A, p0
%             (absolute). Default 1e-3;
% Dsort     : sort model in order of increasing diffusion constants.
%             Default=false.
% display   : Level of output. 0: no output. 1 (default): print convergence
%             message. 2: print convergence every iteration.
%
% 2016-06-28 : researched convergence problems, and found one that was due
% to one state becoming completely unoccopied, which in turn induces NaNs
% in the transition matrix, and from there to the whole model.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMhmm.MLEconverge, variational diffusive HMM, postion estimator
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


% default parameter values
lnLrelTol=1e-8;
parTol=1e-3;
Nwarmup=5;
maxIter=5000;
showConv_lnL=false;
showExit=true;
sortModel=false;

% parameter interpretations
nv=1;
while(nv <= length(varargin))
   pname=varargin{nv};
   if(~ischar(pname))
       error(['optinal arguments must be name/value pairs.'])
   end
   pval=varargin{nv+1};
   nv=nv+2;
    
   if(strcmp(pname,'Nwarmup'))
      Nwarmup=pval;
   elseif(strcmp(pname,'maxIter'))
      maxIter=pval;      
   elseif(strcmp(pname,'lnLrelTol'))
      lnLrelTol=pval;      
   elseif(strcmp(pname,'parTol'))
      parTol=pval;      
   elseif(strcmp(pname,'Dsort'))
      sortModel=pval;    
   elseif(strcmp(pname,'display'))
      n=pval;
      switch n
          case 0
              showConv_lnL=false;
              showExit=false;
          case 1
              showConv_lnL=false;
              showExit=true;
          case 2
              showConv_lnL=true;
              showExit=true;
          otherwise
              error(['Did not understand display ' int2str(n)])
      end
              
   else
       error(['Unrecognized option ' pname ])
   end   
end

% construct convergence report
EMexit=struct;
EMexit.stopcondition='maxiter';
% convergence iterations
lnL0=-inf;
EMtimer=tic;

converged_lnL=0;
converged_par=false;
dParam=inf;
dlnLrel=inf;
for r=1:(Nwarmup+maxIter)
    if(sortModel)
        % sort in order of increasing diffusion constant
        [~,ind]=sort(W.P.lambda);
        if(prod(ind==sort(ind))==0) % the resort the model
            W.P.lambda=W.P.lambda(ind);
            W.P.p0=W.P.p0(ind);
            W.P.A=W.P.A(ind,ind);
            W.S.pst=W.S.pst(:,ind);
            W.S.wA=W.S.wA(ind,ind);
        end
    end

    % iterate
    W=EMhmm.hiddenStateUpdate(W,dat);
    dlnLrel=(W.lnL-lnL0)/abs(W.lnL);
    lnL0=W.lnL;
    W=EMhmm.diffusionPathUpdate(W,dat);
    
    if(r>Nwarmup)
        W=EMhmm.MLEparameterUpdate(W,dat);
        if(exist('lam0','var'))
            dLam=max(abs(W.P.lambda-lam0)./W.P.lambda);
            dA=max(abs(W.P.A(:)-A0(:)));
            dp0=max(abs(W.P.p0-p00));
            dParam=max([ dLam dA dp0]);
            if(r>(Nwarmup+2))
                if(dParam<parTol && ~converged_par)
                    converged_par=true;
                    EMexit.stopcondition='parTol';
                end
            end
        end
        % parameter convergence
        lam0=W.P.lambda;
        A0=W.P.A;
        p00=W.P.p0;
    end
    
    if(showConv_lnL)
        disp(['it ' int2str(r) ', dlnL = ' num2str(dlnLrel,4) ', dPar = ' num2str(dParam,4) ])
    end
    if(r>(Nwarmup+2) && abs(dlnLrel)<lnLrelTol && converged_lnL<4)
        converged_lnL=converged_lnL+1;
        EMexit.stopcondition='lnLrelTol';
    else
        converged_lnL=0;
    end
    if(converged_lnL>=4 && converged_par)
        break
    end
end
EMexit.time=toc(EMtimer);

% add convergence report to model struct
EMexit.numiter=r;
EMexit.dlnLrel=dlnLrel;
W.EMexit=EMexit;
if(showExit)
    disp(W.EMexit)
end

%% path estimates
if(nargout>=2) % compute sequence of most likely states
    [W,sMaxP]=EMhmm.hiddenStateUpdate(W,dat);
end
if(nargout>=3) % compute Viterbi path, with a small offset to avoid infinities
    [W,sMaxP,sVit]=EMhmm.hiddenStateUpdate(W,dat);
end
