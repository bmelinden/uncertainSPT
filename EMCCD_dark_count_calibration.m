function [gain,sigRead,offset,dc,F]=EMCCD_dark_count_calibration(IM,gain0,sigRead0,dc0,doPlot)
% [gain,sigRead,offset,dc]=EMCCD_dark_count_calibration(IM,gain0,sigRead0,dc0,doPlot)
% Calibration from dark counts movie, using the high gain approximation
% with uniform dark current and Gaussian readout noise. The model assumes
% that the dark current is low, meaning dc<<1. In many cases, such images
% can be aquired by closing the camera shutter.
%
% Statistical model: each pixel is modeled as a constant offset, iid
% Gaussian noise with std sigRead, and (with proibability dc) a
% contribution from a dark current electron, which is assumed to follow the
% same EMCCD statistics as normal photoelectrons. I.e., 
% IM(i,j,t) = offset(i,j)+sigRead*n(i,j,t)+s(i,j,t)*m(i,j,t), where
% offset  = pixel-dependent camera offset
% sigRead = std of electronic readous noise 
% n(i,j,t)= iid unit Gaussian random variables
% s(i,j,t)= iid indicator variables (=0,1) for photon present (s=1, with
% probability dc) or abent (s=0, probability 1-dc).
%
% Input:
% IM    : dark movie to calibrate from. A large movie converges slower, but
%         one needs quite many counts in order to determine gain from the
%         high tail of the movie.
% gain0     : initial guess for gain (AD-corrected)
% sigRead0  : initial guess for readout noise standard deviation
%             Good initial guesses inproves runtime and convergence.
% dc0       : initial guess for dark current probability.
%             Good initial guesses inproves runtime and convergence.
% doPlot    : optional figure number for plotting calibration fit (default
%             0, no plot) 
% output 
% gain, sigRead : gain ad readout noise parameters
% offset    : offset map (double)
% dc        : average dark photon count per pixel and frame, should be <<1
%             for the calibration model (at most 1 photon in each pixel) to
%             be a valid approximation. 
% F         : the converged log-likelihood
% Commandline output monitors the convergence (i=iteration number):
% sigRead(i) gain(i) dc(i) dPrel (F(i)-F(i-1))/|F(i)|,
% where dPrel is the maximum relative change of sigRead, gain, mean(dc).
% Iterations terminate when dPrel < 1e-5.
% ML 2016-04-21

% parameters
if(~exist('doPlot','var'))
    doPlot=0;
end
%% initial guess
qq=dc0;
ss=sigRead0;
gg=gain0;
oo=median(IM,3);
oot=repmat(oo,1,1,size(IM,3));

%% plot distributions
if(doPlot>0)
    figure(doPlot)
    clf
    pause(0.1)
    dIM=IM-oot;
    [a,cBins]=hist(dIM(:),min(dIM(:)):1:max(dIM(:)));
    db=mean(diff(cBins));
    a=a/sum(a)/db;
    a=a(2:end-1);
    cBins=cBins(2:end-1);
    
    aMin=min(a(a>0));
    figure(doPlot)
    clf
    subplot(2,1,1)
    hold on
    imHist=stairs(cBins-db/2,a,'k','linew',2);
    set(gca,'yscale','log')
    
    % MLE fit parameters
    lnp0=-log(ss)-0.5*log(2*pi)-0.5*(cBins/ss).^2+log(1-qq);
    lnp1=-log(2)-log(gg)+ss^2/2/gg^2-cBins/gg...
        +log(1+erf(cBins/ss/sqrt(2)-ss/gg/sqrt(2)))+log(qq);
    pFit=exp(lnp0)+exp(lnp1);
    imFit=plot(cBins,pFit,'r-','linew',1);
    
    legend('histogram','MLE fit')
    xlabel('count-offset')
    ylabel('log(pdf)')
    box on
    axis([-5*ss max(dIM(:)) aMin/2 2*max(a(a>0))])
    
    tstr=title('');
    
    subplot(2,1,2)
    imO=imagesc(oo);
    axis equal
    axis tight
    set(gca,'clim',[-10+min(oo(:)) 10+max(oo(:))])
    colormap gray
    title('offset')
    colorbar
    pause(0.1)
end
%% EM iterations

optmin=optimset(@fminunc);
optmin=optimset(optmin,'display','off');%,'algorithm','quasi-newton');
optsol=optimset(@fsolve);
optsol=optimset(optsol,'display','off');%,'algorithm','quasi-newton');
format compact
disp('starting iterations')
dPrel=1;
F=-inf;
disp('       rms |      gain |      d.c. |    dP/|P| | dlnL/|lnL| |         ')
while(dPrel>1e-5)
    
    % start parameters
    sgq0=[ss gg qq];
    
    % responsibilities update
    lnp0=-log(ss)-0.5*log(2*pi)-0.5*((IM-oot)/ss).^2+log(1-qq);
    lnp1=-log(2)-log(gg)+ss^2/2/gg^2-(IM-oot)/gg...
        +log(1+erf( (IM-oot)/ss/sqrt(2)-ss/gg/sqrt(2)))+log(qq);
    Q=1./(1+exp(lnp0-lnp1));

    Fold=F;
    lnpMax=max(lnp0,lnp1);
    F=sum(lnpMax(:)+log(exp(lnp0(:)-lnpMax(:))+exp(lnp1(:)-lnpMax(:))));    
    dFrel=(F-Fold)/abs(F);

    
    % parameter update
    qq=mean(Q(:));
        
    % offset: fixed-point iterations
    dOmax=1;
    it=0;
    maxIter=10;
    while(dOmax>1e-8)
        % attempt fixed-point iterations
        % summary statistics
        T0=sum(1-Q,3);
        T1=sum(Q,3);
        ck0=sum((1-Q).*IM,3)./T0;
        
        oot=repmat(oo,1,1,size(IM,3));
        % iterative update
        oo1=ck0+T1./T0*ss^2/gg...
            -ss./T0*sqrt(2/pi).*sum(...
            Q.*exp(-((IM-oot)/ss/sqrt(2)-ss/gg/sqrt(2)).^2)...
            ./(1+erf((IM-oot)/ss/sqrt(2)-ss/gg/sqrt(2))),3);
        % 
        dOmax=max(abs(oo1(:)-oo(:)));
        oo=oo1;
        it=it+1;
        if( it>maxIter )
            disp(['max(do_k) = ' num2str(dOmax)  ' not converged.'])
            break
        end
    end
    oot=repmat(oo,1,1,size(IM,3));
    dIM=IM-oot;
    
    % gain&noise: fixed-point iterations
    dsdg=1;
    it=0;
    maxIter=10;
    while(dsdg>1e-8)
        T1    = sum(Q(:));
        T0    = sum(1-Q(:));
        cmo2_0= sum( (1-Q(:)).*(dIM(:)).^2)/T0;
        cmo_1 = sum( Q(:).*(dIM(:)))/T1;
        
        W1    = sum( Q(:).*dIM(:).*sqrt(2/pi).*...
            exp(-0.5*(dIM(:)/ss-ss/gg).^2)./(...
            1+erf((dIM(:)/ss-ss/gg)/sqrt(2))));
        W2    = sum( Q(:).*sqrt(2/pi).*...
            exp(-0.5*(dIM(:)/ss-ss/gg).^2)./(...
            1+erf((dIM(:)/ss-ss/gg)/sqrt(2))));
        
        sg0=[ss gg];
        sgFun=@(sg)([-1+cmo2_0/sg(1)^2+(sg(1)/sg(2))^2*T1/T0-W1/T0/sg(1)-sg(1)/sg(2)*W2/T0;
            -1-(sg(1)/sg(2))^2+cmo_1/sg(2)+sg(1)/sg(2)*W2/T1  ]);
        sg1=fsolve(sgFun,sg0,optsol);
        dsdg=max(abs( (sg1-sg0)./sg1));
        %disp(num2str([sg1 dgds]))
        ss=sg1(1);
        gg=sg1(2);
        it=it+1;
        if( it>maxIter )
            disp(['max(d(s,g)/|s,g|) = ' num2str(dsdg)  ' not converged.'])
            break
        end
    end
    sgq=[ss gg qq];
    dPrel=max(abs((sgq-sgq0)./sgq));
    disp(num2str([ss gg qq dPrel dFrel]))
    if(doPlot>0)
        db=mean(diff(cBins));
        [a,cBins]=hist(dIM(:),[cBins(1)-db cBins cBins(end)+db]);
        a=a/sum(a)/db;
        a=a(2:end-1);
        cBins=cBins(2:end-1);
        set(imHist,'ydata',a);
        
        % MLE fit parameters
        lnp0=-log(ss)-0.5*log(2*pi)-0.5*(cBins/ss).^2+log(1-qq);
        lnp1=-log(2)-log(gg)+ss^2/2/gg^2-cBins/gg...
            +log(1+erf(cBins/ss/sqrt(2)-ss/gg/sqrt(2)))+log(qq);
        pFit=exp(lnp0)+exp(lnp1);
        set(imFit,'ydata',pFit);
        set(imO,'cdata',oo)
        
        set(tstr,'String',sprintf('gain=%.1f, rms=%.1f, dark curr.=%.1e',gg,ss,qq))
        pause(0.1)
    end
end
gain=gg;
sigRead=ss;
offset=oo;
dc=qq;
   
     



