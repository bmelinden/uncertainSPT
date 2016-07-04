function W=MLEmodelsearch(dat,R,tau,N,varargin)
% W=MLEmodelsearch(dat,R,tau,N,varargin)
% Perform a model search for an optimal MLE fit of an EMhmm model to data.
%
% Output:
% W     : optimal HMM model
% if asked for, the fields W.sMaxP and/or W.sVit are added.
% sMaxP : trajectory of most likely states 
% sVit  : most likely trajectory of hidden states (from Viterbi algorithm)
%
% Input :
% dat   : EMhmm data struct, e.g., from EMhmm.preprocess
% R,tau : blur factors for dat
% N     : number of hidden states
% additional arguments in the form of 'name', value pairs
%
% restarts  : number of MLE search cycles to run. Minimum 1, Default 12.
% wA        : pseudocount matrix for guessing a transition matrix A : 
%             A(k,:) ~ Dir(wA(k,:)). Default: wA=ones(N,N)+2*eye(N).
% Ddt       : range of diffusion constant*dt-values to use for initial
%             guesses. Default [500 50000]
% viterbi   : true/false. Compute viterbi path estimate of hidden states
%             for the optimal model, save in W.sVit. Default: false.
% sMaxP     : true/false. Compute the path of most likely states (not the
%             same as the Viterbi path), save in W.sMaxP. Default: false.
% lnLrelTol : relative convergence criteria for (lnL(n)-lnL(n-1))/|lnL(n)|.
%             Default 1e-8;
% parTol    : convergence criteria for parameters lambda (relative), A, p0
%             (absolute). Default 1e-3;
% maxIter   : maximum number of iterations (past Nwarmup). Default 5000.
% scalePar  : [scaleMin scaleSteps]. scaleMin is the smallest scaling value
%             to use for rms/blur parameters. scaleSteps is the number of
%             scaling steps. Default : [1e-3 5];
%
% Name-value pairs not recognized a passed on to EMhmm.MLEconverge at
% appropriate times, e.g.,
%
% Nwarmup   : number of initial iterations where model parameters are kept
%             constant in order to 'burn in' the states and hidden path.
%             Default 5.
% Dsort     : sort model in order of increasing diffusion constants.
%             Default=false.
% display   : Level of output. 0: no output. 1 (default): print convergence
%             message. 2: print convergence every iteration.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMhmm.MLEmodelsearch
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


%% parameter interpretations
% default parameter values
lnLrelTol=1e-8;
parTol=1e-3;
wA=ones(N,N)+2*eye(N);
logDdt=log([5e2 5e4]);
restarts=12;
doViterbi=false;
doSmaxP=false;
scaleMin=1e-2;
scaleSteps=7;
maxIter=5000;

showConv_lnL=false;
showExit=true;

% default MLEconverge options
MLEargs={'Dsort',true,'display',0,'lnLrelTol',lnLrelTol,'parTol',parTol};

% parse arguments
nv=1;
while(nv <= length(varargin))
   pname=varargin{nv};
   if(~ischar(pname))
       error(['optinal arguments must be name/value pairs.'])
   end
   pval=varargin{nv+1};
   nv=nv+2;
   
   if(strcmp(pname,'restarts'))
       restarts=pval;
   elseif(strcmp(pname,'wA'))
       wA=pval;
       if(size(wA,1) ~= N )
           error('Pseudocount matrix size incompatible with number os states.')
       end
   elseif(strcmp(pname,'Ddt'))
       logDdt(1)=log(pval(1));
       logDdt(2)=log(pval(2));
       logDdt=sort(logDdt);
   elseif(strcmp(pname,'Viterbi') || strcmp(pname,'viterbi') )
       doViterbi=pval;
   elseif(strcmp(pname,'SmaxP'))
       doSmaxP=pval;
   elseif(strcmp(pname,'lnLrelTol'))
       lnLrelTol=pval;
       MLEargs{end+1}='lnLrelTol';
       MLEargs{end+1}=lnLrelTol;
   elseif(strcmp(pname,'parTol'))
       parTol=pval;
       MLEargs{end+1}='parTol';
       MLEargs{end+1}=parTol;
   elseif(strcmp(pname,'maxIter'))
       maxIter=pval;
       MLEargs{end+1}='maxIter';
       MLEargs{end+1}=maxIter;
   elseif(strcmp(pname,'scalePar'))
       scaleMin=pval(1);
       scaleSteps=pval(2);
   else
       disp(['Passing unrecognized option ' pname ' to MLEconverge.'])
       MLEargs{end+1}=pname;
       MLEargs{end+1}=pval;
   end   
end
%% model search
% no error data
dat0=dat;
dat0.v=dat0.v*scaleMin^2; % almost not errors accounted for

% construct convergence report
ModSearchExit=struct;
ModSearchExit.method='';

%% loop over restarts, MLE and no blur-err init
Wmle=cell(1,restarts);
Wnb=cell(1,restarts);
%disp('starting MLE pafor loop')
for n=1:restarts
    % construct initial model
    A_init=dirrnd(wA);
    p0_init=A_init^1000;
    p0_init=p0_init(1,:);
    Ddt_init=exp(sort(logDdt(1)+rand(1,2*N-1)*diff(logDdt)));
    Ddt_init=Ddt_init(1:2:end); % some spread between diffusion constants    
    % simple MLE convergence    
    tconv=tic;
    Wmle{n}=EMhmm.init_P_dat(tau,R,Ddt_init,A_init,p0_init,dat);   
    disp(['converging MLE model ' int2str(n)])% ' of ' int2str(restarts) ])
    W0=Wmle{n};
    Wmle{n}=EMhmm.MLEconverge(Wmle{n},dat,MLEargs{:});
    disp(['converged MLE model ' int2str(n) ', ' int2str(size(A_init,1)) ' states, ' num2str(toc(tconv)/60)  ' min.'])% ' of ' int2str(restarts) ])
    pause(0.1)
    Wmle{n}.modSearch='mle rnd init';
    % initial guess based on no-blur-fit
    tic
    Wnb{n}=EMhmm.init_P_dat(0,0,Ddt_init,A_init,p0_init,dat0);    
    Wnb{n}=EMhmm.MLEconverge(Wnb{n},dat0,MLEargs{:});
end
% best MLE fit
MLElnL=zeros(1,restarts);
for n=1:restarts
    MLElnL(n)=Wmle{n}.lnL;
end
[~,b]=max(MLElnL);
W=Wmle{b};
% best no-blur-no err fit
WNB=struct;
WNB.lnL=-inf;
for n=1:restarts
    if(Wnb{n}.lnL>WNB.lnL)
        WNB=Wnb{n};
    end
end
WSC=WNB; % save for later use in the scaling approach
%% converge from no-blu-no-err initialization
WNB.tau=tau;
WNB.R=R;
iter=0;
while(iter<maxIter)
    lnL0=WNB.Y.MeanLnQ;
    WNB=EMhmm.diffusionPathUpdate(WNB,dat);
    dlnL=abs((WNB.Y.MeanLnQ-lnL0)/WNB.lnL);
    
    lam0=WNB.P.lambda;
    A0=WNB.P.A;
    p00=WNB.P.p0;
    WNB=EMhmm.MLEparameterUpdate(WNB,dat);
    dLam=max(abs(WNB.P.lambda-lam0)./WNB.P.lambda);
    dA=max(abs(WNB.P.A(:)-A0(:)));
    dp0=max(abs(WNB.P.p0-p00));
    dParam=max([ dLam dA dp0]);
    
    %disp(num2str([dlnL dParam]))
    iter=iter+1;
    if(iter>3 && dlnL<lnLrelTol && dParam < parTol)
        break
    end
end
WNB=EMhmm.MLEconverge(WNB,dat,MLEargs{:});
WNB.modSearch='no-blur-err init';
if(WNB.lnL>W.lnL)
    W=WNB;
end
%% logarithmic scaling approach
tic
scaleFactor=logspace(log10(scaleMin),0,scaleSteps);
for k=1:numel(scaleFactor)
    F=scaleFactor(k);
    WSC.tau=tau*F;
    WSC.R=R*F;
    dat0=dat;
    dat0.v=dat0.v*F^2;
    WSC=EMhmm.MLEparameterUpdate(WSC,dat0);
    WSC=EMhmm.MLEconverge(WSC,dat0,MLEargs{:},'Nwarmup',0,'display',0);
    %disp(['Done (repl step maxStep) : '  int2str([nrs k numel(scaleFactor)])])
end
WSC.modSearch=['log scale, ' int2str(numel(scaleFactor)) ' steps.'];
if(WSC.lnL>W.lnL)
    W=WSC;
end
%% tidy up the selected model
modSearch.method=W.modSearch;
modSearch.lnL=W.lnL;
modSearch.lnL_nbe=WNB.lnL;
modSearch.lnL_lsc=WSC.lnL;
modSearch.lnL_mle=MLElnL;
W.modSearch=modSearch;
% path estimates
if(doSmaxP||doViterbi)
    [W,sMaxP,sVit]=EMhmm.hiddenStateUpdate(W,dat);
    if(doSmaxP)
        W.sMaxP=sMaxP;
    end
    if(doViterbi)
        W.sVit=sVit;
    end
end
% convergence report
if(showExit)
    disp(W.modSearch)
end
