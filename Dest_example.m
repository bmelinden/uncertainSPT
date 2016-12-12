% this script demonstrates the use ov the various D-estimators in the +Dest
% folder. For details on the various files, see Dest/README
clear
%% generate some data
p0=1;
A=1;
pE=0;
D=1e6;
dt=10e-3;
tE=10e-3;

Ddt=D*dt;
lambda_tru=2*D*dt;
R  =tE/dt/6;
tau=tE/dt/2;
beta=tau*(1-tau)-R;

locErrG=[10 2];    % localization error std s~Gam, with <s>=20, std(s)=5 

T=100;
numTrj=100;
dim=1;

[x,v,~,y]=EMhmm.diffusiveHMM_blur_detach(p0,A,pE,Ddt,tE/dt,locErrG,dim,T,numTrj);

if(~iscell(x))
    x={x};
    v={v};
    y={y};
end

% add explicit errors estimates
xv=cell(size(x));
for k=1:length(xv)
    xv{k}=[x{k} sqrt(v{k})];
end

%dat=ML1_preprocess(xv,dim);
dat=Dest.preprocess_mixed_columns(xv,1:dim,(dim+1):(2*dim),[],false);
                                
[logL,trj]=Dest.logLlambda(dat,lambda_tru,tau,R);

%% plot some stuff
figure(1)
clf
hold on
box on

hh=errorbar(1:T,dat.x(dat.one(1):dat.end(1)),xv{1}(:,2),'k');

hh(3)=plot(1:T+1,trj.mu(dat.Yone(1):dat.Yend(1)),'r','linew',1);
plot(1:T+1,trj.mu(dat.Yone(1):dat.Yend(1))+sqrt(trj.CovDiag0(dat.Yone(1):dat.Yend(1))),'r--')
plot(1:T+1,trj.mu(dat.Yone(1):dat.Yend(1))-sqrt(trj.CovDiag0(dat.Yone(1):dat.Yend(1))),'r--')
hh(2)=plot(y{1},'-b.','linew',2);

legend(hh,'data\pmstd.','true positions','est. positions\pm std')
xlabel('time')
ylabel('position')
%% estimate lambda using ML estimator

f=@(L)(-Dest.logLlambda(dat,L,tau,R));

[lambda_ML,fMin]=fminsearch(f,lambda_tru);
DML=lambda_ML/2/dt;

ll=lambda_tru*logspace(-0.3,0.3,30);
logL=ll;
for k=1:length(ll)
    logL(k)=-f(ll(k));
end

figure(2)
clf
hold on
plot(ll/lambda_tru,logL,'-k.')
plot(lambda_ML/lambda_tru,-fMin,'r*')
grid on

xlabel(' D / D_{true}')
ylabel('ln L')
set(gca,'xscale','log')
