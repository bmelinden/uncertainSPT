function W=MLEparameterUpdate(W,dat)
% one round of EM iterations in a diffusive HMM, including handling missing
% data, and adding conjgate prior fields
% P0.wA, P0.w0, P0.n, P0.c if the model W has not prior field P0.

beta=W.tau*(1-W.tau)-W.R;
tau=W.tau;

%% parameter update
W.P.A=rowNormalize(W.S.wA);
W.P.p0=rowNormalize(sum(W.S.pst(W.one,:),1));

% second attempt at diffusion constants
MLEopt = optimoptions('fminunc','GradObj','off','TolX',1e-14,...
    'MaxIter',1e5,'TolFun',1e-14,'Algorithm','quasi-newton',...
    'MaxFunEvals',1e10,'DerivativeCheck','off','Display','none');

indS=1:sum(dat.T+1);
indS(W.end)=0;
indS=indS(indS>0);
dy2=sum((W.Y.mu(indS+1,:)-W.Y.mu(indS,:)).^2 ...
    +W.Y.sig0(indS,:)+W.Y.sig0(indS+1,:)-2*W.Y.sig1(indS,:),2);
dxy2=(dat.x(indS,:)-W.Y.mu(indS,:)*(1-tau)-W.Y.mu(indS+1,:)*tau).^2 ...
        +(1-tau)^2*W.Y.sig0(indS,:)...
        +tau^2*W.Y.sig0(indS+1,:)...
        +2*tau*(1-tau)*W.Y.sig1(indS,:);
dxy2(~isfinite(dat.v(indS,:)))=0; %%% ad-hoc warning, use pO instead
for j=1:W.N
    pt=W.S.pst(indS,j);
    nj=W.dim/2*sum(pt);    
    cj=0.5*sum(pt.*dy2);
    
    Fj=@(L)(-nj*log(L)-cj/L-0.5*sum(pt.*sum(log1p(beta*L./dat.v(indS,:))+dxy2./(dat.v(indS,:)+beta*L),2)));
    
    mLogLj=@(logLam)(-Fj(exp(logLam)));
    W.P.lambda(j)=exp(fminunc(mLogLj,log(W.P.lambda(j)),MLEopt));
end



