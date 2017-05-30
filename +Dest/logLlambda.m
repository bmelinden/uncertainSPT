function [logL,trj]=logLlambda(dat,lambda,tau,R)
%  [logL,trj]=Dest.logLlambda(dat,lambda,tau,R)
% log likelihood of diffusion constant for trj with varying localization
% errors
%
% input:
% dat       : trajectory data struct from Dest.preprocess_mixed_columns
% lambda    : step vaiance lambda=2*D*dt
% tau       : exposure average time. 
% R         : Berglund's blur factor. 
%             For illumination with illumination fraction
%             tExposure/tSample, we have 
%             tau = tExposure/tSample/2 for constant
%             R = tExposure/tSample/6
%
% output:
% logL      : log likelihood of diffusion constant, with hidden positions
%             integrated out
% trj       : struct with information about the hidden trajectory.
% trj.mu    : expected positions mu(t) = <y(t)>
% trj.one   : sub-trajectory start indices
% trj.end   : sub-trajectory end indices
%             Sub-trajectory i stretches from
%             trj.mu(trj.one(i):trj.end(i)), and has length T(i)+1, if
%             there were T(i) measured positions in trj i.
% trj.CovDiag0 : diagonal covariance matrix entries, 
%                trj.CovDiag0(t) = <(y(t)-<y(t)>)^2>
% trj.CovDiag1 : first off-diagonal covariance matrix entries, 
%                trj.CovDiag1(t) = <(y(t)-<y(t)>)*(y(t+1)-<y(t+1)>)>
% Note that more long-range correlations also exist, but are not computed.
%             
% Note: if you want to look under the hood, note that these functions may
% use partly different index conventions than those in +EMhmm, and expect 
% precision as STANDARD DEVIATIONS, while EMhmm work with variances.

%% start of actual code
%% preprocess the data and insert default aggregation
beta=tau*(1-tau)-R;
if(beta<0)
    error('ML1_logLlambda: inconsistent blur parameters ( tau*(1-tau)-R<0 =')
end
%% compute likelihood
%% assemble the Lambda matrix
% construct \nu
Ttot=sum(dat.T+1);   % total number of steps in hidden trajectory
nu=zeros(Ttot,1);
Ldiag=zeros(Ttot,3); % diagonals in the Lambda matrix
BL=beta*lambda;
for nt=1:length(dat.Yone)
    MU1=dat.Yone(nt);
    MUe=dat.Yend(nt)-1; % 1:T for hidden trajectory, which has T+1 points
    X1=dat.one(nt);     % 1,T for the measured trajectory, which has T points
    XT=dat.end(nt);
    T=1+XT-X1; % length of current trajectory
    vTot_t=BL+dat.v(X1:XT,1); % total variance in this trj
    
    
    nu(MU1,1)      =dat.x(X1,1)/vTot_t(1)*(1-tau);
    nu(MUe+1,1)    =dat.x(XT,1)/vTot_t(end)*tau;
    nu(MU1+1:MUe,1)=dat.x(X1+1:XT,1)./vTot_t(2:end)*(1-tau)...
                   +dat.x(X1:XT-1,1)./vTot_t(1:end-1)*tau;
    
    % diagonal elements
    Ldiag(MU1,2)      =1/lambda+(1-tau)^2/vTot_t(1);
    Ldiag(MU1+1:MUe,2)=2/lambda+(1-tau)^2./vTot_t(2:end)+tau^2./vTot_t(1:end-1);
    Ldiag(MUe+1,2)    =1/lambda+tau^2/vTot_t(end);
    % off-diagonal elements
    Ldiag(MU1:MUe,1)    =tau*(1-tau)./vTot_t-1/lambda;
    Ldiag(MU1+1:MUe+1,3)=tau*(1-tau)./vTot_t-1/lambda;
end
% compute \Sigma(t,t), \Sigma(t,t+1)
%Lambda=spdiags(Ldiag,-1:1,Ttot,Ttot);
% invert Lambda block by block
% Sigma1=inv(Lambda);
% Sigma=spalloc(Ttot,Ttot,sum((dat.Yend-dat.Yone+1).^2));


%%%% do it in one go with a custom algorithm? Yes!
[CovDiag0,CovDiag1,logDetLambda]=...
    triSym_triInv_rescale_trjWise(Ldiag(:,2),Ldiag(:,1),dat.Yone,dat.Yend,length(dat.Yone));
%    triSym_d1Inv_trjWise(Ldiag(:,2),Ldiag(:,1),dat.Yone,dat.Yend,length(dat.Yone));

% one giant linear system?
Lambda=spdiags(Ldiag,-1:1,size(Ldiag,1),size(Ldiag,1));
mu=Lambda\nu;

dmunux2=0;
for nt=1:length(dat.Yone)
    MU1=dat.Yone(nt);
    MUe=dat.Yend(nt)-1; % 1:T for hidden trajectory, which has T+1 points
    X1=dat.one(nt);     % 1,T for the measured trajectory, which has T points
    XT=dat.end(nt);
    
    % doing the quadratic contribution while minimizing differences btw
    % large numbers
    stdTot=sqrt(BL+dat.v(X1:XT,1));
    munu=sqrt(mu(MU1:MUe+1).*nu(MU1:MUe+1));
    dmunux2=dmunux2+0.5*sum( (munu(1:end-1)-dat.x(X1:XT)./stdTot )...
                           .*(munu(1:end-1)+dat.x(X1:XT)./stdTot ))...        
                            +0.5*munu(end).^2;    
end

% assemble total likelihood
logL=-sum(dat.T)/2*log(lambda)-0.5*sum(log(dat.v(:,1)+BL))-0.5*logDetLambda+dmunux2;

% collect averages for hidden trajectory
trj.mu=mu;
trj.CovDiag0=CovDiag0;
trj.CovDiag1=CovDiag1;


