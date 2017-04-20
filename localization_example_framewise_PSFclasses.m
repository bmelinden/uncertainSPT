% An example of using the EMCCDfit code to track a single particle
% and estimate uncertainties with and without prior. This uses a
% frame-by-frame analysis, which is perhaps more efficient when analyzing
% large image sets. Also, a set of PSF classes that incorporates both PSF
% model and prior distributions are used. 
% ML 2016-08-31

clear
addpath(pwd)

% parameters
EMgain=90;
sigmaRead=20;
ROIwidth=9;
Nquad=5;

% MLE 
PsgMLE=PSF.SymGauss_MLE('initialGuess',[0 0 log([1 100 1.5])],'priorParameters',[]);
PagMLE=PSF.AsymGauss_angle_MLE('initialGuess',[0 0 log([1 100 1.25 1.75]) 0],'priorParameters',[]);

% MAP
PsgMAP=PSF.SymGaussS0_logNormBN_expS0('initialGuess',[0 0 log([1 100 0.1])],...
            'priorParameters',[0 log(30) log(200) inf 1],'NA',1.4,'lambda',639/80); % wavelength 639 nm, px size 80 nm
PagMAP=PSF.AsymGaussS0_BNlnN_expS0('initialGuess',[0 0 log([1 100 0.75 0.5]) 0],...
            'priorParameters',[0 2 2 2 1],'NA',1.4,'lambda',639/80);

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
        [frameCoord,frameCov,framePar]=EMCCDfit.refineSingleFrame(...
            frameCoord0,currentFluoFrame,fluoOffset,ROIwidth,Nquad,logLobj,PsgMAP);
        fitTime(1)=fitTime(1)+toc(t0);
        
        coord_sg_MAP(ind,:)=frameCoord; % save frame number
        param_sg_MAP(ind)=framePar;
        xyCov_sg_MAP(ind,:)=frameCov;
        
        % symmetric Gaussian, MLE fit
        t0=tic;
         [frameCoord,frameCov,framePar]=EMCCDfit.refineSingleFrame(...
            frameCoord0,currentFluoFrame,fluoOffset,ROIwidth,Nquad,logLobj,PsgMLE);
        fitTime(2)=fitTime(2)+toc(t0);

        coord_sg_MLE(ind,:)=frameCoord; % save frame number
        param_sg_MLE(ind)=framePar;
        xyCov_sg_MLE(ind,:)=frameCov;       
        
        % asymmetric Gaussian, MAP fit
        t0=tic;
        [frameCoord,frameCov,framePar]=EMCCDfit.refineSingleFrame(...
            frameCoord0,currentFluoFrame,fluoOffset,ROIwidth,Nquad,logLobj,PagMAP);
        fitTime(3)=fitTime(3)+toc(t0);

        coord_ag_MAP(ind,:)=frameCoord; % save frame number
        param_ag_MAP(ind)=framePar;
        xyCov_ag_MAP(ind,:)=frameCov;

        % asymmetric Gaussian, MLE fit
        t0=tic;
        [frameCoord,frameCov,framePar]=EMCCDfit.refineSingleFrame(...
            frameCoord0,currentFluoFrame,fluoOffset,ROIwidth,Nquad,logLobj,PagMLE);
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
subplot(2,1,1)
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

subplot(2,1,2)
hold on
plot(sqrt(xyCov_ag_MAP(:,1)),'-b')
plot(sqrt(xyCov_ag_MLE(:,1)),'-r')
xlabel('frame')
ylabel('Laplace uncertainty [px]')
legend('MAP Laplace','MLE Laplace')
title('asymmetric Gauss uncertainty')
box on


ind=find(sqrt(xyCov_ag_MLE(:,1))<=2*a);

disp('Error and uncertainties (excluding points with est. Laplace uncert. >2 px) :')
disp('------------------------------')
disp('symmetric Gaussian PSF model:')
ind=find(sqrt(xyCov_sg_MAP(:,1))<=2*a);
disp(['MAP retained: ' num2str(numel(ind)/T,4)])
disp(['MAP RMS err.: ' num2str(sqrt(mean((coord_sg_MAP(ind,1)-emTrj_px(ind,1)).^2)),4) ' px'])
disp(['MAP Laplace : ' num2str(sqrt(mean(xyCov_sg_MAP(ind,1))),4) ' px'])
disp(['MAP CRLB    : ' num2str(sqrt(mean(MAPxCRLB(ind))),4) ' px'])
disp('-')
ind=find(sqrt(xyCov_sg_MLE(:,1))<=2*a);
disp(['MLE retained: ' num2str(numel(ind)/T,4)])
disp(['MLE RMS err.: ' num2str(sqrt(mean((coord_sg_MLE(ind,1)-emTrj_px(ind,1)).^2)),4) ' px'])
disp(['MLE Laplace : ' num2str(sqrt(mean(xyCov_sg_MLE(ind,1))),4) ' px'])
disp(['MLE CRLB    : ' num2str(sqrt(mean(MLExCRLB(ind))),4) ' px'])
disp('------------------------------')
disp('asymmetric Gaussian PSF model:')
ind=find(sqrt(xyCov_ag_MAP(:,1))<=2*a);
disp(['MAP retained: ' num2str(numel(ind)/T,4)])
disp(['MAP RMS err.: ' num2str(sqrt(mean((coord_ag_MAP(ind,1)-emTrj_px(ind,1)).^2)),4) ' px'])
disp(['MAP Laplace : ' num2str(sqrt(mean(xyCov_ag_MAP(ind,1))),4) ' px'])
disp('-')
ind=find(sqrt(xyCov_ag_MLE(:,1))<=2*a);
disp(['MLE retained: ' num2str(numel(ind)/T,4)])
disp(['MLE RMS err.: ' num2str(sqrt(mean((coord_ag_MLE(ind,1)-emTrj_px(ind,1)).^2)),4) ' px'])
disp(['MLE Laplace : ' num2str(sqrt(mean(xyCov_ag_MLE(ind,1))),4) ' px'])
disp('------------------------------')
