% a simple example of ifferent ways to run the HMM algorithm
% ML 2016-07-04
addpath(genpath(pwd)); % add this folder and subfolders to matlab path

% parameters
p0=ones(3,1)/3;
A=rowNormalize(diag(2*[15 30 40])+ones(3,3));
D=[0.1 1 3]*1e6;  % nm^2/s
dt=20e-3;            % s
tE=dt;              % exposure time
pM=0.1;
locErrG=[20 10];    % localization error std s~Gam, with <s>=20, std(s)=10 
dim=2;
T=[150 100];

% blur coefficients
tau=tE/dt/2;
R=tE/dt/6;

% generate synthetic data
[x,v,s,y]=EMhmm.diffusiveHMM_blur_detach(p0,A,0,D*dt,tE/dt,locErrG,dim,T,numel(T),pM);

% preprocess
dat=EMhmm.preprocess(x,v,2,s);

%% analysis 1: a single initial guess
% initialize an HMM model with the true parameters
W1=EMhmm.init_P_dat(tau,R,D*dt,A,p0,dat);
W1=EMhmm.MLEconverge(W1,dat,'Dsort',true,'Nwarmup',10);

%% analysis 2: search for global optimum 
W2=EMhmm.MLEmodelsearch(dat,R,tau,W1.N,'Dsort',true,'Ddt',[0.1 5]*1e6*dt);

%% analysis 3: cheat using the known ground truth
W3=W1;
W3.P.lambda=D*dt;
W3.P.p0=ones(1,W3.N)/W3.N;
W3.P.A=A;
W3.S.pst=W3.S.pst*0;
for k=1:numel(dat.misc)
    if(dat.misc(k)>0)
        W3.S.pst(k,dat.misc(k))=1;
    end
end

% start with on diffusion path update since parameters and hidden states
% are completely correct to start with
W3=EMhmm.diffusionPathUpdate(W3,dat);
% then run standard convergence
W3=EMhmm.MLEconverge(W3,dat);

disp('log likelihood (highest wins) :')
fprintf('init w true parameters: lnL = %.2f \n',W1.lnL)
fprintf('model search          : lnL = %.2f \n',W2.lnL)
fprintf('ground truth init     : lnL = %.2f \n',W3.lnL)


%% plot mean diffusion constant vs time

Dt=dat.misc;
Dt(Dt>0)=D(dat.misc(Dt>0))/1e6;
D1=W1.S.pst*W1.P.lambda'/2/dt/1e6; % [um^2/s]
D2=W2.S.pst*W2.P.lambda'/2/dt/1e6; % [um^2/s]
D3=W3.S.pst*W3.P.lambda'/2/dt/1e6; % [um^2/s]

figure(1)
clf
hold on

% do not plot points btw trajectories
Dt(Dt==0)=nan;
D1(D1==0)=nan;
D2(D2==0)=nan;
D3(D3==0)=nan;

plot(Dt,'k','linew',2)
plot(D1,'b')
plot(D2,'r')
plot(D3,'g--')

legend('truth','parameter init','model search','ground truth init')
xlabel('frame')
ylabel('D [\mum^2/2]')

%% perform a simple bootstrap analysis on the 'model search' model
Nbs=100;
Pmle=EMhmm.parameterEstimate(W2,dt);
BS=EMhmm.parameterBootstrap(W2,dat,Nbs,dt,true);
EMhmm.displayParameterBootstrap(Pmle,BS);

