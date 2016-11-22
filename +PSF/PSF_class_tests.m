% test script for PSF classes
clear
meta.abstractDetails('PSF.PSFModel')
meta.abstractDetails('PSF.AsymGauss_angle')
meta.abstractDetails('PSF.SymGauss')

% create an object and test functions

P=PSF.AsymGauss_angle_logNormB('initialGuess',[3 3 log([1 200 1.75 2.25]) 0],'priorParameters',[0 1])
[E,dEdp]=P.psfDensity(0,0,P.initialGuess);
[y,dy]=P.logPrior(P.initialGuess);
par0 = P.convertToOutStruct(P.initialGuess, struct)

P=PSF.AsymGauss_angle_MLE('initialGuess',[3 3 log([1 200 1.75 2.25]) 0],'priorParameters',[])
[E,dEdp]=P.psfDensity(0,0,P.initialGuess);
[y,dy]=P.logPrior(P.initialGuess);
par0 = P.convertToOutStruct(P.initialGuess, struct)

P=PSF.SymGauss_MLE('initialGuess',[3 3 log([1 200 1.75])],'priorParameters',[])
[E,dEdp]=P.psfDensity(0,0,P.initialGuess);
[y,dy]=P.logPrior(P.initialGuess);
par0 = P.convertToOutStruct(P.initialGuess, struct)

P=PSF.SymGauss_logNormB('initialGuess',[3 3 log([1 200 1.75])],'priorParameters',[0 1])
[E,dEdp]=P.psfDensity(0,0,P.initialGuess);
[y,dy]=P.logPrior(P.initialGuess);
par0 = P.convertToOutStruct(P.initialGuess, struct)

P=PSF.AsymGaussS0_MLE('initialGuess',[0 0 log([1 200 0.75 0.5]) 0],'priorParameters',[],'NA',1.4,'lambda',639/80)
[E,dEdp]=P.psfDensity(0,0,P.initialGuess);
[y,dy]=P.logPrior(P.initialGuess);
[par0,f0] = P.convertToOutStruct(P.initialGuess, struct)

P=PSF.AsymGaussS0_BNlnN_dSexp_N('initialGuess',[0 0 log([1 200 0.75 0.5]) 0],'priorParameters',[0 2 2 2 1 log(1.5)],'NA',1.4,'lambda',639/80)
[E,dEdp]=P.psfDensity(0,0,P.initialGuess);
[y,dy]=P.logPrior(P.initialGuess);
[par0,f0] = P.convertToOutStruct(P.initialGuess, struct)

% log-prior derivatives



% derivatives of intensity model
k=6;
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


subplot(2,1,2)
plot(pk,dlnP0dk,'linew',2)
hold on
plot((pk(1:end-1)+pk(2:end))/2,diff(lnP0)./diff(pk),'r--','linew',2)




