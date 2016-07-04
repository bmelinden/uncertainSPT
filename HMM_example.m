% a simple example of how to run the HMM algorithm

addpath(genpath(pwd)); % add this folder and subfolders to matlab path

% parameters
p0=ones(3,1)/3;
A=rowNormalize(diag([30 20 15])+ones(3,3));
dt=5e-3;            % s
D=[0.1 0.5 4]*1e6;  % nm^2/s
tE=dt;              % exposure time
pM=0.1;
locErrG=[20 10];    % localization error std s~Gam, with <s>=20, std(s)=10 
dim=2;
T=[150 100];

% blur coefficients
tau=tE/dt/2;
R=tE/dt/6;

% generate synthetic data
[x,v,s,y]=EMhmm.diffusiveHMM_blur_detach(p0,A,0,D*dt,tE,locErrG,dim,T,numel(T),pM);

% preprocess
dat=EMhmm.preprocess(x,v,2);

% analysis 1: a single initial guess
% initialize an HMM model
W=EMhmm.init_P_dat(tau,R,D*dt,A,p0,dat);
W=EMhmm.MLEconverge(W,dat,'Dsort',true);



% plot results


