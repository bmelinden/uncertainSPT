% test script for PSF classes, to 
% - illustrate how to create objects
% - verify that analytical partical derivatives of PSF intensaity models
%   and log_prior densities are implemented correctly

% PSF master class
meta.abstractDetails('PSF.PSFModel') 

% derived PSF models, no priors implemented
meta.abstractDetails('PSF.SymGauss')
meta.abstractDetails('PSF.AsymGauss_angle')
meta.abstractDetails('PSF.SymGaussS0')
meta.abstractDetails('PSF.AsymGaussS0')

classNumber=7 % which PSF class to test and verify

% create PSF objects from some PSF & prior classes
switch classNumber
    case 1
        P=PSF.SymGauss_MLE('initialGuess',[3 3 log([1 200 1.75])],'priorParameters',[])
        [E,dEdp]=P.psfDensity(0,0,P.initialGuess);
        [y,dy]=P.logPrior(P.initialGuess);
        par0 = P.convertToOutStruct(P.initialGuess, struct)
    case 2
        P=PSF.SymGauss_logNormB('initialGuess',[3 3 log([1 200 1.75])],'priorParameters',[0 log(20)])
        [E,dEdp]=P.psfDensity(0,0,P.initialGuess);
        [y,dy]=P.logPrior(P.initialGuess);
        par0 = P.convertToOutStruct(P.initialGuess, struct)
    case 3
        P=PSF.SymGauss_logNormBN('initialGuess',[3 3 log([1 200 1.75])],...
            'priorParameters',[0 log(20) log(200) log(20)])
        [E,dEdp]=P.psfDensity(0,0,P.initialGuess);
        [y,dy]=P.logPrior(P.initialGuess);
        par0 = P.convertToOutStruct(P.initialGuess, struct)
    case 4
        P=PSF.SymGauss_logNormBNS('initialGuess',[3 3 log([1 200 1.75])],...
            'priorParameters',[0 log(20) log(200) log(20) log(0.21*639/1.4/80) log(20)])
        [E,dEdp]=P.psfDensity(0,0,P.initialGuess);
        [y,dy]=P.logPrior(P.initialGuess);
        par0 = P.convertToOutStruct(P.initialGuess, struct)
    case 5
        P=PSF.AsymGauss_angle_MLE('initialGuess',[3 3 log([1 200 1.75 2.25]) 0],'priorParameters',[])
        [E,dEdp]=P.psfDensity(0,0,P.initialGuess);
        [y,dy]=P.logPrior(P.initialGuess);
        par0 = P.convertToOutStruct(P.initialGuess, struct)
    case 6
        P=PSF.SymGaussS0_MLE('initialGuess',[0 0 log([1 200 0.1])],...
            'priorParameters',[],'NA',1.4,'lambda',639/80) % wavelength 639 nm, px size 80 nm
        [E,dEdp]=P.psfDensity(0,0,P.initialGuess);
        [y,dy]=P.logPrior(P.initialGuess);
        [par0,f0] = P.convertToOutStruct(P.initialGuess, struct)
    case 7
        P=PSF.SymGaussS0_expS0('initialGuess',[0 0 log([1 200 0.1])],...
            'priorParameters',[1],'NA',1.4,'lambda',639/80) % wavelength 639 nm, px size 80 nm
        [E,dEdp]=P.psfDensity(0,0,P.initialGuess);
        [y,dy]=P.logPrior(P.initialGuess);
        [par0,f0] = P.convertToOutStruct(P.initialGuess, struct)
    case 8
        P=PSF.SymGaussS0_logNormBN_expS0('initialGuess',[0 0 log([1 200 0.1])],...
            'priorParameters',[0 log(30) log(200) log(30) 1],'NA',1.4,'lambda',639/80) % wavelength 639 nm, px size 80 nm
        [E,dEdp]=P.psfDensity(0,0,P.initialGuess);
        [y,dy]=P.logPrior(P.initialGuess);
        [par0,f0] = P.convertToOutStruct(P.initialGuess, struct)
    case 9
        P=PSF.AsymGaussS0_MLE('initialGuess',[0 0 log([1 200 0.75 0.5]) 0],...
            'priorParameters',[],'NA',1.4,'lambda',639/80)
        [E,dEdp]=P.psfDensity(0,0,P.initialGuess);
        [y,dy]=P.logPrior(P.initialGuess);
        [par0,f0] = P.convertToOutStruct(P.initialGuess, struct)
    case 10
        P=PSF.AsymGaussS0_BNlnN_expS0('initialGuess',[0 0 log([1 200 0.75 0.5]) 0],...
            'priorParameters',[0 2 2 2 1],'NA',1.4,'lambda',639/80)
        [E,dEdp]=P.psfDensity(0,0,P.initialGuess);
        [y,dy]=P.logPrior(P.initialGuess);
        [par0,f0] = P.convertToOutStruct(P.initialGuess, struct)
end

%% derivatives of intensity model
% check that analytical and numerical partial derivatives match
for k=1:numel(P.initialGuess)
    p0=P.initialGuess;
    pk=p0(k)+linspace(-5,5,5000);
    E=0*pk;
    dEdpk=E;
    xx=0.14;
    yy=0.3;
    
    for m=1:numel(pk)
        p1=p0;
        p1(k)=pk(m);
        [E(m),dEdp]= P.psfDensity(xx, yy, p1);
        dEdpk(m)=dEdp(k);
        [lnP0(m),dlnP0] = P.logPrior(p1);
        dlnP0dk(m)=dlnP0(k);
    end
    
    figure(1)
    clf
    
    subplot(2,1,1)
    plot(pk,dEdpk,'linew',2)
    hold on
    plot((pk(1:end-1)+pk(2:end))/2,diff(E)./diff(pk),'r--','linew',2)
    title(['intensity derivative, parameter no. ' int2str(k)])
    legend('analytical','numerical','location','best')
    
    subplot(2,1,2)
    plot(pk,dlnP0dk,'linew',2)
    hold on
    plot((pk(1:end-1)+pk(2:end))/2,diff(lnP0)./diff(pk),'r--','linew',2)
    title(['log-prior derivative, parameter no. ' int2str(k)])
    legend('analytical','numerical','location','best')

    pause(0.5)
    
end



