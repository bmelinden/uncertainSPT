% a script demonstrating the EMCCD calibration function on simulated data.
% Edit below to change simulation paramters.

% offset with some striped noise
offset=100+ones(120,1)*round(5*rand(1,500))+round(1.5*randn(120,500));
sigRead=15;
gain=60;
q=0.01; % expected photons/pixel
T=100;

disp('creating test data...')
% not quite uniform intentisy
qMap=q*ones(size(offset));
Nph=poissrnd(repmat(qMap,1,1,T));
IM=round(gamrnd(Nph,gain)+sigRead*randn(size(Nph)));
IM=IM+repmat(offset,1,1,size(IM,3));
disp('running EMCCDcal.dark_count_calibration')
[gainCal,sigReadCal,offsetCal,dcCal]=EMCCD_dark_count_calibration(IM,gain,sigRead,q,1);

% compre input and output
disp(['EM gain      : calibrated ' num2str(gainCal,3) ', true value ' num2str(gain,3)])
disp(['readout noise: calibrated ' num2str(sigReadCal,4) ', true value ' num2str(sigRead,4)])
disp(['dark current : calibrated ' num2str(dcCal  ,3) ', true value ' num2str(q   ,3)])
disp(['offset       : calibration error mean=' num2str(mean(offsetCal(:)-offset(:)),3) ', std=' ...
            num2str(std(offsetCal(:)-offset(:)),3)])
        



