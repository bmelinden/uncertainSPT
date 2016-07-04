% An example of using the EMCCDfit code to track a single particle
% and estimate uncertainties with and without prior.
% ML 2016-07-04

clear
addpath(genpath(pwd))

% parameters
EMgain=90;
sigmaRead=20;
ROIwidth=9;

% construct a localization prior
lnBGmu=log(0.8); % weak background prior since background is time-dependent
lnBGstd=2;
Smu0=log(1.5);   % same skew-normal distribution as used in main text
Sstd0=1;
Salpha=5;

% prior for symmetric Gaussian PSF model
lnPrior_lnBNS_sym=@(lnBNS)(-0.5*(lnBNS(1)-lnBGmu)^2/lnBGstd^2 ...
            +EMCCDfit.skewGauss_logPdf(lnBNS(3),Smu0,Sstd0,Salpha));
fminOpt = optimoptions('fminunc','TolX',1e-6,'MaxIter',500,...
    'TolFun',1e-6,'MaxFunEvals',10000,'Display','notify','algorithm','quasi-newton' );

% indata: this is a simulated movie, and thus we use the known true
% positions instead of spot detection to create an initial guess
MV=EMCCDfit.ML_loadStack2('trickyEvent438491736029240_01.tif');
fluoOffset=double(imread('fake_offset_50x50.tif'));
R=load('trickyEvent438491736029240.mat');

% extract true coordinates
emTrj_nm=R.emissionAverage{1};
nm2px=1/R.opt.camera.pixLength;
emTrj_px=[emTrj_nm(:,1:3)*nm2px emTrj_nm(:,4:end)]; % the initial 'guess'
clear R

% construct an MLEobject
if(0)
    cMax=max(MV(:))-min(fluoOffset(:));
    Emax=cMax/EMgain*5;
    logLobj=EMCCDfit.logL_EMCCD_lookup(EMgain,sigmaRead,cMax,Emax);
    save logL.mat logLobj cMax Emax EMgain sigmaRead
else
    load logL.mat
end

%% localize all positions

T=size(emTrj_nm,1);
MAPcoord=zeros(T,3);
MAPbNS  =zeros(T,3);
MAPxyCov=zeros(T,4);
MLEcoord=zeros(T,3);
MLEbNS  =zeros(T,3);
MLExyCov=zeros(T,4);
fprintf('spot refinement:             ')
for t=1:size(emTrj_px,1)
    tic
    frame=emTrj_px(t,4);    
    dotCoord=emTrj_px(t,[1 2 4]);
    currentFluoFrame=MV(:,:,frame);

    % construct ROI for localization
    x0   =dotCoord(1);
    y0   =dotCoord(2);
    
    [riSpot,ciSpot,x0Spot,y0Spot]=EMCCDfit.ROItransform(x0,y0,ROIwidth,ROIwidth,size(currentFluoFrame));
    spotROI=double(currentFluoFrame(riSpot,ciSpot)); % movie raw data
    bgROI=double(fluoOffset(riSpot,ciSpot));  % offset background to be subtracted
    fitData=(spotROI-bgROI);                  % daata to fit to

    % MAP fit
    MAPobj_sym=EMCCDfit.logL_psf(logLobj,@EMCCDfit.psf_diff_symgauss,fitData,3);
    MAPfun_sym=@(p)(-MAPobj_sym.lnL(p)-lnPrior_lnBNS_sym(p(3:5)));  % -log(likelihood)-log(prior)
    lnpInit=[x0-x0Spot y0-y0Spot log(10) log(300) log(2)];      % initial guess in the ROI
    [lnpMAPdot,~,exitFlag,~,~,hessianMAP] = fminunc(MAPfun_sym,lnpInit,fminOpt);
    covMAP=inv(hessianMAP); % numerical covariance matrix from 
    hessRcond=rcond(hessianMAP);

    % transform back to global coordinates
    xMAP=x0Spot+lnpMAPdot(1);
    yMAP=y0Spot+lnpMAPdot(2);
    % other fit parameters
    bg=exp(lnpMAPdot(3));
    N =exp(lnpMAPdot(4));
    S =exp(lnpMAPdot(5));
    tFit=toc;
    
    MAPcoord(t,1:2)=[xMAP yMAP];
    MAPcoord(t,3)=dotCoord(3); % save frame number
    MAPbNS(t,:)  =[bg N S];
    MAPxyCov(t,:)=[covMAP(1,1) covMAP(1,2) covMAP(2,2) hessRcond];
    
    % MLE fit
    MLEobj_sym=EMCCDfit.logL_psf(logLobj,@EMCCDfit.psf_diff_symgauss,fitData,3);
    MLEfun_sym=@(p)(-MLEobj_sym.lnL(p));  % -log(likelihood)
    lnpInit=[x0-x0Spot y0-y0Spot log(10) log(300) log(2)];      % initial guess in the ROI
    [lnpMLEdot,~,exitFlag,~,~,hessianMLE] = fminunc(MLEfun_sym,lnpInit,fminOpt);
    covMLE=inv(hessianMLE); % numerical covariance matrix from 
    hessRcond=rcond(hessianMLE);

    % transform back to global coordinates
    xMLE=x0Spot+lnpMLEdot(1);
    yMLE=y0Spot+lnpMLEdot(2);
    % other fit parameters
    bg=exp(lnpMLEdot(3));
    N =exp(lnpMLEdot(4));
    S =exp(lnpMLEdot(5));
    tFit=toc;
    
    MLEcoord(t,1:2)=[xMLE yMLE];
    MLEcoord(t,3)=dotCoord(3); % save frame number
    MLEbNS(t,:)  =[bg N S];
    MLExyCov(t,:)=[covMLE(1,1) covMLE(1,2) covMLE(2,2) hessRcond];
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b %4d %4d ',t,T);
end
fprintf('\n')

% scale back to nm
MAPcoord(:,1:2)=MAPcoord(:,1:2)/nm2px;
MAPbNS(:,3)=MAPbNS(:,3)/nm2px;
MAPxyCov(:,1:3)=MAPxyCov(:,1:3)/nm2px^2;

MLEcoord(:,1:2)=MLEcoord(:,1:2)/nm2px;
MLEbNS(:,3)=MLEbNS(:,3)/nm2px;
MLExyCov(:,1:3)=MLExyCov(:,1:3)/nm2px^2;

%% compute CRLB
a=1/nm2px;

% MAP fit
s2a=MAPbNS(:,3).^2+a^2/12;
N  =MAPbNS(:,2);
b  =MAPbNS(:,1);
tt =2*pi*s2a.*b./N/a^2;

MAPxCRLB= 2*s2a./N.*(1+4*tt+sqrt(2*tt./(1+4*tt)));

% MLE fit
s2a=MLEbNS(:,3).^2+a^2/12;
N  =MLEbNS(:,2);
b  =MLEbNS(:,1);
tt =2*pi*s2a.*b./N/a^2;

MLExCRLB= 2*s2a./N.*(1+4*tt+sqrt(2*tt./(1+4*tt)));

%% compare estimated RMS uncertainties
figure(1)
clf
hold on


plot(sqrt(MAPxyCov(:,1)),'-b')
plot(sqrt(MLExyCov(:,1)),'-r')
plot(sqrt(MAPxCRLB(:,1)),'b.')
plot(sqrt(MLExCRLB(:,1)),'r.')
xlabel('frame')
ylabel('Laplace uncertainty [nm]')
legend('MAP Laplace','MLE Laplace','MAP CRLB','MLE CRLB')
box on

ind=find(sqrt(MAPxyCov(:,1))<=2*a);


disp('Error and uncertainties (excluding points with RMS uncert. >2 px) :')
disp(['MAP RMS err.: ' num2str(sqrt(mean((MAPcoord(ind,1)-emTrj_nm(ind,1)).^2)),4) ' nm'])
disp(['MAP Laplace : ' num2str(sqrt(mean(MAPxyCov(ind,1))),4) ' nm'])
disp(['MLE Laplace : ' num2str(sqrt(mean(MLExyCov(ind,1))),4) ' nm'])
disp(['MAP CRLB    : ' num2str(sqrt(mean(MAPxCRLB(ind,1))),4) ' nm'])
disp(['MLE CRLB    : ' num2str(sqrt(mean(MLExCRLB(ind,1))),4) ' nm'])
    




