function [dotCoord,dotCov,dotParam,rawParam]=...
    refineSingleFrame(dotCoord0,fluoFrame,fluoOffset,ROIwidth,Nquad,logLobj,psfObj,fitParam0,saveError)
% [dotCoord,dotCov,dotParam,rawParam]=...
%    refineSingleFrame(dotCoord0,fluoFrame,fluoOffset,ROIwidth,Nquad,...
%    logLobj,psfObj,fitParam0,saveError)
%
% attempt MAP fit of the suggested dots in dotCoord0 to a single
% fluorescent image fluoFrame, with mean offset fluoOffset.
%
% dotCoord0=[x y ...]. Suggested starting points for MLE fits. Columns > 2
%                      are just passed on unchanged to dotCoordMLE.
% fluoFrame : image to fit spots to.
% fluoOffset: EMCCD camera offset image
% ROIwidth  : side length (pixels) of the ROI square to use for fit
% Nquad     : number of quadrature point per dimension for estimating total
%             pixel intensity
% logLobj   : logL_EMCCD_lookup object (this is how EMgain and readoutNoise
%             are included in the fit).
% psfObj    : psf object with prior and initial guess, subclass of PSF.PSFmodel 
% fitParam0 : optional complete initial guess array. Replaces the initial
%             guess that is normally derived from dotCoord0 and the psfObj.
%             Useful for multi-step optimization, when the output of one
%             optimization step is used as initial guess for the next.
% saveError : Optional. If true, save workspace for points that somehow
%             fail to converge to disk. Useful for debugging. 
%             Default: false.
%
% dotCoord  : same as dotCoord0, but with first two columns updated to x,y
%             positions from the fit.
% dotCov    : position covariances, dotCov(t,:)= [varX covXY varY],
%             from a Laplace approximation.
% dotParam  : a parameter struct array, with fields specified by the
%             psf2param function+additional fields
%             dr = distance between refined position and initial guess
%             logL = MAP fit log(likelihood)
%             detH = |determinant of logL Hessian|
% rawParam  : raw fit parameter array, as returned by fminunc; this is the
%             same format as that of fitParam0
%
% For spots where the optimization fails, dotCov(t,:)=[inf inf inf]
%
% Martin Lindén, bmelinden@gmail.com, 2016-09-26
externalInitialGuess=false;
if(exist('fitParam0','var') && ~isempty(fitParam0))    
    externalInitialGuess=true;
end
if(~exist('saveError','var') || isempty(saveError))    
    saveError=false;
end

dotCoord=dotCoord0;
dotCov=zeros(size(dotCoord,1),3);
p0=psfObj.convertToOutStruct(psfObj.initialGuess,struct);
p0.dr=NaN;             % add refinement displacement
p0.logL=NaN;
p0.logLaic=NaN;
p0.logLbic=NaN;
dotParam(1:size(dotCoord,1))=p0;
clear p0

fminOpt = optimoptions(EMCCDfit.pointOptimization.fminuncDefault,...
    'TolX',1e-6,'MaxIter',500,'TolFun',1e-6,'MaxFunEvals',10000,'Display','notify' );

for r=1:size(dotCoord0,1)
    tic
    % construct ROI for localization
    x0   =dotCoord0(r,1);
    y0   =dotCoord0(r,2);    
    [riSpot,ciSpot,x0Spot,y0Spot]=EMCCDfit.ROItransform(x0,y0,ROIwidth,ROIwidth,size(fluoFrame));
    spotROI=double(fluoFrame(riSpot,ciSpot)); % movie raw data
    bgROI=double(fluoOffset(riSpot,ciSpot));  % offset background to be subtracted
    fitData=(spotROI-bgROI);
    
    % MAP fit
    MAPobj=EMCCDfit.pointOptimization(logLobj,psfObj,fitData,Nquad); % likelihood function, data, and numerical quadrature
    lnpInit=psfObj.initialGuess;
    lnpInit(1:2)=[x0-x0Spot y0-y0Spot];      % initial guess in the ROI

    if(externalInitialGuess)
        lnpInit=fitParam0(r,:);
    end
    if(r==1)
        rawParam=zeros(size(dotCoord0,1),size(lnpInit,2));
    end
    rawParam(r,:)=lnpInit;
    try
        [lnpMAPdot,logLMAP,exitFlag,~,~,hessianMAP] = ...
            fminunc(@MAPobj.minus_lnL,lnpInit,fminOpt);
        covMAP=inv(hessianMAP); % numerical covariance matrix
        hessRcond=rcond(hessianMAP);
        detTrace=[det(covMAP(1:2,1:2)) trace(covMAP(1:2,1:2))];
        rawParam(r,:)=lnpMAPdot;
        logDetH=sum(log(eig(hessianMAP)));        
    catch me
        %warning('MAP optimization error! Writing NaN results.')
        %disp(me)
        % if something goes wrong already here, we discard
        % the point
        exitFlag=-inf;
        logLMAP=inf;
        lnpMAPdot=nan(1,length(lnpInit));
        hessianMAP=inf(1,2);
        hessRcond=0;
        covMAP=inf(length(lnpInit),length(lnpInit));
        detTrace=inf(1,2);
        logDetH=0;
    end
    xMAP=x0Spot+lnpMAPdot(1);
    yMAP=y0Spot+lnpMAPdot(2);
    % sanity check
    if( exitFlag<=0 || ~isfinite(hessRcond) || hessRcond<10*eps || sum((detTrace<=0))>0)
        % the there are convergence problems of some kind
        disp('--------------------------------------------------------')
        warning('localization convergence problem or ill-conditioned Hessian:')
        disp(['fminunc exitFlag = ' int2str(exitFlag)])
        disp(['rcond            = ' num2str(hessRcond)])
        disp(['fit parameters   = ' num2str(lnpMAPdot,4)])
        disp('parameter struct :')
        psfObj.convertToOutStruct(lnpMAPdot)
        %disp(['eig (covXY)      = ' num2str(eig(covMAP(1:2,1:2))',4)])
        disp('PSF model:')
        psfObj
        disp('--------------------------------------------------------')
        if(saveError)
            if(exist('errLogs','dir')==7 || mkdir('errLogs'))
                errFile=sprintf('errLogs%serrLog_refineSingleFrame_%09d.mat',filesep,round(1e9*rand));
            else
                errFile=sprintf('errLog_refineSingleFrame_%09d.mat',round(1e9*rand));
            end
            save(errFile)
        end
        if(0)% a visual debug block to visualize data and failed fits
            figure(13)
            clf
            subplot(2,1,1)
            hold on
            imagesc(fluoFrame)
            plot(xMAP,yMAP,'ro')
            plot(lnpInit(1)+x0Spot,lnpInit(2)+y0Spot,'gs')
            
            axis equal, axis tight
            cr=get(gca,'clim');
            
            subplot(2,2,3)
            hold on
            imagesc(fitData+bgROI)
            plot(lnpMAPdot(1),lnpMAPdot(2),'ro')
            plot(lnpInit(1),lnpInit(2),'gs')
            axis equal, axis tight
            set(gca,'clim',cr);
            title('ROI data')
            
            subplot(2,2,4)
            hold on
            Emodel=MAPobj.psfModel(lnpMAPdot);
            imagesc(Emodel*logLobj.EMgain+bgROI)
            plot(lnpMAPdot(1),lnpMAPdot(2),'ro')
            plot(lnpInit(1),lnpInit(2),'gs')
            axis equal, axis tight
            set(gca,'clim',cr);
            title('ROI fitted model')
            
            colormap gray
        end
        % if the fit failed, mark it as such
        xMAP=nan;
        yMAP=nan;
        covMAP=inf(length(lnpInit),length(lnpInit));
    end
    dotCoord(r,1:2)=[xMAP yMAP];
    dotCov(r,:)=[covMAP(1,1) covMAP(1,2) covMAP(2,2)];
    dotParam(r) = psfObj.convertToOutStruct(lnpMAPdot,dotParam(r));
    
    dr=norm(dotCoord(r,1:2)-dotCoord0(r,1:2));
    dotParam(r).logL=-logLMAP;        
    K=numel(lnpInit);
    dotParam(r).logLaic=dotParam(r).logL-K;
    dotParam(r).logLbic=dotParam(r).logL-1/2*logDetH+K/2*log(2*pi);
    dotParam(r).dr=dr;	% add refinement displacement
    
    if( false && ( dotCov(1,1)<0 || dotCov(2,2)<0 ) )
        % a visual debug block: it seems most fits that encounter this
        % warning are false positives anyway, so the block is not
        % activated.
        disp('--------------------------------------------------------')
        disp(['fminunc exitFlag = ' int2str(exitFlag)])
        disp('hessianMLE : ')
        disp(num2str(hessianMAP,4))
        disp(['rcond = ' num2str(hessRcond)])
        disp('mleData : ')
        disp(num2str(fitData,5))
        disp('fit parameters :')
        disp(num2str(lnpMAPdot,4))
        
        figure(13)
        clf
        subplot(2,1,1)
        hold on
        imagesc(fluoFrame)
        plot(xMAP,yMAP,'ro')
        axis equal, axis tight
        cr=get(gca,'clim');
        
        subplot(2,2,3)
        hold on
        imagesc(fitData)
        plot(lnpMAPdot(1),lnpMAPdot(2),'ro')
        axis equal, axis tight
        set(gca,'clim',cr);
        
        subplot(2,2,4)
        hold on
        Emodel=mapObj.psfModel(lnpMAPdot);
        imagesc(Emodel*logLobj.EMgain)
        plot(lnpMAPdot(1),lnpMAPdot(2),'ro')
        axis equal, axis tight
        set(gca,'clim',cr);
        
        colormap gray
        
        warning('Negative localization error variance encountered')
        disp('--------------------------------------------------------')
    end
    %disp([ int2str([r size(dotCoord0,1)]) ' ' num2str(tFit) ' s'])
end
