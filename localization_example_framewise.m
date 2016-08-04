% An example of using the EMCCDfit code to track a single particle
% and estimate uncertainties with and without prior. This uses a
% frame-by-frame analysis, which is perhaps more efficient when analyzing
% large image sets.
% ML 2016-08-04

clear
addpath(genpath(pwd))

% parameters
EMgain=90;
sigmaRead=20;
ROIwidth=9;
Nquad=5;

% parameters for a localization prior
lnBGmu=log(0.8); % weak background prior since background is time-dependent
lnBGstd=2;
Smu0=log(1.5);   % same skew-normal distribution as used in main text
Sstd0=1;
Salpha=5;

% MLE: flat prior
psfLogNoPrior=@(x)(0);
% prepare for using a symmetric Gaussian PSF model
psfFun_symGauss=@EMCCDfit.psf_diff_symgauss;
psfLogPrior_symGauss=@(xylnBNS)(-0.5*(xylnBNS(3)-lnBGmu)^2/lnBGstd^2 ...
            +EMCCDfit.skewGauss_logPdf(xylnBNS(5),Smu0,Sstd0,Salpha));
psfInit_symGauss=[0 0 log([0.8 100 1.5])];
psf2param_symGauss=@EMCCDfit.psf2param_symgauss;


% prepare for using an asymmetric Gaussian PSF model
psfFun_asymGauss=@EMCCDfit.psf_diff_asymgauss_angle;
% fit parameters : [ muX muY lnB lnN lnS1 lnS2 v]
psfLogPrior_asymGauss=@(xylnBNS)(...
            -0.5*(xylnBNS(3)-lnBGmu)^2/lnBGstd^2 ... % background prior
            +EMCCDfit.skewGauss_logPdf(mean(xylnBNS(5:5)),Smu0,Sstd0,Salpha));% ...
            %+EMCCDfit.skewGauss_logPdf(xylnBNS(6),Smu0,Sstd0,Salpha));
            % ad-hoc prior
psfInit_asymGauss=[0 0 log([0.8 100 1.5 1.6]) 0]; % start with a little bit of asymmetry 
psf2param_asymGauss=@EMCCDfit.psf2param_asymgauss_angle;

% indata: this is a simulated movie, and thus we use the known true
% positions instead of spot detection to create an initial guess
MV=EMCCDfit.ML_loadStack2('trickyEvent438491736029240_01.tif');
fluoOffset=double(imread('fake_offset_50x50.tif'));
R=load('trickyEvent438491736029240.mat');

% extract true coordinates
emTrj_nm=R.emissionAverage{1}; % x y z frame cell  
nm2px=1/R.opt.camera.pixLength;
emTrj_px=[emTrj_nm(:,1:3)*nm2px emTrj_nm(:,4:end)]; % the initial 'guess'
clear R

if(0)% multiply number of points to test indexing and measure time
    emTrj_px=[emTrj_px;emTrj_px];
    emTrj_px=[emTrj_px;emTrj_px];
    emTrj_px=[emTrj_px;emTrj_px];
    emTrj_px=[emTrj_px;emTrj_px];
    emTrj_px=[emTrj_px;emTrj_px];
end

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
coord_sg_MAP=zeros(size(emTrj_px));
coord_sg_MLE=zeros(size(emTrj_px));
coord_ag_MAP=zeros(size(emTrj_px));
coord_ag_MLE=zeros(size(emTrj_px));

clear param_sg_MAP param_sg_MLE param_ag_MAP param_ag_MLE

xyCov_sg_MAP=zeros(T,3);
xyCov_sg_MLE=zeros(T,3);
xyCov_ag_MAP=zeros(T,3);
xyCov_ag_MLE=zeros(T,3);

fitTime=zeros(1,4);
Nfits=0;

%fprintf('spot refinement in frame:             ')
%fprintf('                                    ')
fprintf('        frame ')
fprintf('[ symMAP symMLE asymMAP asymMLE]')
fprintf('\n              ')
fprintf('                                        ')

for frame=1:size(MV,3)
    ind=find(emTrj_px(:,4)==frame);
    if(~isempty(ind));
        currentFluoFrame=double(MV(:,:,frame));
        frameCoord0=emTrj_px(ind,:);
        
        % symmetric Gaussian, MAP fit
        t0=tic;
        [frameCoord,frameCov,framePar]=...
            EMCCDfit.MAP_EMCCD_refineSingleFrame(frameCoord0,currentFluoFrame,fluoOffset,ROIwidth,Nquad,logLobj,...
            psfFun_symGauss,psfLogPrior_symGauss,psfInit_symGauss,psf2param_symGauss);
        fitTime(1)=fitTime(1)+toc(t0);
        
        coord_sg_MAP(ind,:)=frameCoord; % save frame number
        param_sg_MAP(ind)=framePar;
        xyCov_sg_MAP(ind,:)=frameCov;
        
        % symmetric Gaussian, no prior
        t0=tic;
         [frameCoord,frameCov,framePar]=...
            EMCCDfit.MAP_EMCCD_refineSingleFrame(frameCoord0,currentFluoFrame,fluoOffset,ROIwidth,Nquad,logLobj,...
            psfFun_symGauss,psfLogNoPrior,psfInit_symGauss,psf2param_symGauss);
        fitTime(2)=fitTime(2)+toc(t0);

        coord_sg_MLE(ind,:)=frameCoord; % save frame number
        param_sg_MLE(ind)=framePar;
        xyCov_sg_MLE(ind,:)=frameCov;       
        
        % asymmetric Gaussian, MAP fit
        t0=tic;
        [frameCoord,frameCov,framePar]=...
            EMCCDfit.MAP_EMCCD_refineSingleFrame(frameCoord0,currentFluoFrame,fluoOffset,ROIwidth,Nquad,logLobj,...
            psfFun_asymGauss,psfLogPrior_asymGauss,psfInit_asymGauss,psf2param_asymGauss);
        fitTime(3)=fitTime(3)+toc(t0);

        coord_ag_MAP(ind,:)=frameCoord; % save frame number
        param_ag_MAP(ind)=framePar;
        xyCov_ag_MAP(ind,:)=frameCov;

        % asymmetric Gaussian, MLE fit
        t0=tic;
        [frameCoord,frameCov,framePar]=...
            EMCCDfit.MAP_EMCCD_refineSingleFrame(frameCoord0,currentFluoFrame,fluoOffset,ROIwidth,Nquad,logLobj,...
            psfFun_asymGauss,psfLogNoPrior,psfInit_asymGauss,psf2param_asymGauss);
        fitTime(4)=fitTime(4)+toc(t0);

        coord_ag_MLE(ind,:)=frameCoord; % save frame number
        param_ag_MLE(ind)=framePar;
        xyCov_ag_MLE(ind,:)=frameCov;
        
        Nfits=Nfits+numel(ind);
    end
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b')
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    fprintf('%04d of %04d ',frame,size(MV,3));
    fprintf('[%7.4f%7.4f%8.4f%8.4f] s/spot.',fitTime/Nfits);
end
fprintf('\n')
fitTime=fitTime/Nfits; % average time to fit a single dot
fprintf('\n')

%% compute CRLB for symmetric gauss
a=1; % doing math in pixel units

% MAP fit
s2a=[param_sg_MAP(:).std].^2+a^2/12;
N  =[param_sg_MAP(:).amplitude];
b  =[param_sg_MAP(:).background];
tt =2*pi*s2a.*b./N/a^2;
MAPxCRLB= 2*s2a./N.*(1+4*tt+sqrt(2*tt./(1+4*tt)));

% MLE fit
s2a=[param_sg_MLE(:).std].^2+a^2/12;
N  =[param_sg_MLE(:).amplitude];
b  =[param_sg_MLE(:).background];
tt =2*pi*s2a.*b./N/a^2;
MLExCRLB= 2*s2a./N.*(1+4*tt+sqrt(2*tt./(1+4*tt)));
%% compare estimated RMS uncertainties
figure(1)
clf
hold on

plot(sqrt(xyCov_sg_MAP(:,1)),'-b')
plot(sqrt(xyCov_sg_MLE(:,1)),'-r')
plot(sqrt(MAPxCRLB),'b.')
plot(sqrt(MLExCRLB),'r.')
xlabel('frame')
ylabel('Laplace uncertainty [px]')
legend('MAP Laplace','MLE Laplace','MAP CRLB','MLE CRLB')
title('symmetric Gauss uncertainty')
box on

disp('Error and uncertainties (excluding points with RMS uncert. >3 px) :')
disp('------------------------------')
disp('symmetric Gaussian PSF model:')
ind=find((sqrt(xyCov_sg_MAP(:,1))<=3*a).*(sqrt(xyCov_sg_MLE(:,1))<=3*a));
disp(['MAP RMS err.: ' num2str(sqrt(mean((coord_sg_MAP(ind,1)-emTrj_px(ind,1)).^2)),4) ' px'])
disp(['MAP Laplace : ' num2str(sqrt(mean(xyCov_sg_MAP(ind,1))),4) ' px'])
disp(['MAP CRLB    : ' num2str(sqrt(mean(MAPxCRLB(ind))),4) ' px'])
disp(['MLE RMS err.: ' num2str(sqrt(mean((coord_sg_MLE(ind,1)-emTrj_px(ind,1)).^2)),4) ' px'])
disp(['MLE Laplace : ' num2str(sqrt(mean(xyCov_sg_MLE(ind,1))),4) ' px'])
disp(['MLE CRLB    : ' num2str(sqrt(mean(MLExCRLB(ind))),4) ' px'])
disp('------------------------------')
disp('asymmetric Gaussian PSF model:')
ind=find((sqrt(xyCov_ag_MAP(:,1))<=3*a).*(sqrt(xyCov_ag_MLE(:,1))<=3*a));
disp(['MAP RMS err.: ' num2str(sqrt(mean((coord_ag_MAP(ind,1)-emTrj_px(ind,1)).^2)),4) ' px'])
disp(['MAP Laplace : ' num2str(sqrt(mean(xyCov_ag_MAP(ind,1))),4) ' px'])
disp(['MLE RMS err.: ' num2str(sqrt(mean((coord_ag_MAP(ind,1)-emTrj_px(ind,1)).^2)),4) ' px'])
disp(['MLE Laplace : ' num2str(sqrt(mean(xyCov_ag_MLE(ind,1))),4) ' px'])
%disp(['MAP CRLB    : ' num2str(sqrt(mean(MAPxCRLB(ind))),4) ' px'])
%disp(['MLE CRLB    : ' num2str(sqrt(mean(MLExCRLB(ind))),4) ' px'])
disp('------------------------------')
