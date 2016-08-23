% test script for PSF classes

meta.abstractDetails('SymGauss_logNormB')
meta.abstractDetails('SymGauss_MLE')
meta.abstractDetails('AsymGauss_angle_logNormB')
meta.abstractDetails('AsymGauss_angle_MLE')

% create an object and test functions

P1=SymGauss_logNormB('initialGuess',[3 3 log([1 200 2])],'priorParameters',[0 1])
[E,dEdp]=P1.psfDensity(0,0,P1.initialGuess);
[y,dy]=P1.logPrior(P1.initialGuess);
par1 = P1.convertToOutStruct(P1.initialGuess, struct)

P1=SymGauss_MLE('initialGuess',[3 3 log([1 200 2])])
[E,dEdp]=P1.psfDensity(0,0,P1.initialGuess);
[y,dy]=P1.logPrior(P1.initialGuess);
par1 = P1.convertToOutStruct(P1.initialGuess, struct)

P2=AsymGauss_angle_logNormB('initialGuess',[3 3 log([1 200 1.75 2.25]) 0],'priorParameters',[0 1])
[E,dEdp]=P2.psfDensity(0,0,P2.initialGuess);
[y,dy]=P2.logPrior(P2.initialGuess);
par2 = P2.convertToOutStruct(P2.initialGuess, struct)

P2=AsymGauss_angle_MLE('initialGuess',[3 3 log([1 200 1.75 2.25]) 0])
[E,dEdp]=P2.psfDensity(0,0,P2.initialGuess);
[y,dy]=P2.logPrior(P2.initialGuess);
par2 = P2.convertToOutStruct(P2.initialGuess, struct)
