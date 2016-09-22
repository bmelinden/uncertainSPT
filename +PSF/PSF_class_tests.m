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
