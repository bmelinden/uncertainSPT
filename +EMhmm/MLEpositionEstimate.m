function xEst=MLEpositionEstimate(W,dat)
% xEst=MLEpositionEstimate(W,dat)
% Estimate actual positions (posterior mean position) based on input data
% and converge HMM model.
%
% ML 2016-06-07

beta=W.tau*(1-W.tau)-W.R;
tau=W.tau;

bvl=zeros(size(W.S.pst,1),1); % beta*lambda_t/v_t
ym=zeros(size(W.S.pst,1),1);  % y_t*(1-tau)+y_{t+1}*tau
xEst=zeros(size(W.S.pst,1),W.dim); 

% v 2
for d=1:W.dim
   ym(1:end-1,1)=(1-tau)*W.Y.mu(1:end-1,d)+tau*W.Y.mu(2:end,d);
   bvl=(beta./dat.v(:,d))*W.P.lambda; 
   xEst(:,d)=ym.*sum(W.S.pst./(1+bvl),2)...
            +dat.x(:,d).*sum(W.S.pst.*bvl./(1+bvl),2);
end
xEst(W.end,:)=0;               % make sure there is no position estimate for t=T+1

% variance of position estimate
if(0)
vxEst=zeros(size(xEst));
for d=1:W.dim
   
    % conditional variance
    vzt=(dat.v(:,d)*beta*W.P.lambda)./(...
        ones(size(dat.v(:,d)))*beta*W.P.lambda...
        +dat.v(:,1)*ones(size(W.P.lambda)));
    
    mvzt=mean(vzt.*W.S.pst,2);

    % contribution from hidden position uncertainty and states
    varYm=W.Y.sig0(1:end-1,d)*(1-tau)^2+W.Y.sig0(2:end,d)*tau^2 ...
    +2*tau*(1-tau)*W.Y.sig1(1:end-1,d);


    % contribution from hidden states
    

end

end