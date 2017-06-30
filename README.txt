---------------------------------------------------------------------
uncertainSPT
---------------------------------------------------------------------
uncertainSPT contains two basic software tools to extract and use
localization uncertainty for single particle tracking and localization
microscopy.
1) EMCCDfit, localization algorithms to estimate particle positions
and localization uncertainty from images aquired by an EMCCD camera,
1b) EMCCD_dark_count_calibration.m, a calibration routine that extract
EMCCD parameters from very low light images (<<1 photons/pixel), and
2) EMhmm, a variational EM algorithm that performs maximum likelihood
inference in a diffusive hidden Markov model, where both motion blur
and localization uncertainty are included in the model, through an
extension of the Berglund model for single-state diffusion [1].

The code runs on matlab, with some inner loops implemented in
C/C++. Binaries for 64 bit linux and max OS are included.

If you use this code, please cite our work [2].

1. Berglund, A.J. (2010). Statistics of camera-based single-particle
tracking. Phys. Rev. E 82, 11917. doi: 10.1103/PhysRevE.82.011917

2. Lindén, M., Ćurić, V., Amselem, E., and Elf, J. (2017). Pointwise
error estimates in localization microscopy. Nat Commun 8, 15115.
doi: 10.1038/ncomms15115
======================================================================
Copyright (C) 2016 Martin Lindén and Johan Elf

E-mail: bmelinden@gmail.com, johan.elf@icm.uu.se
======================================================================
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or any
later version.  This program is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
=====================================================================
This product includes code developed by others, see copyright in
individual files:
dirrnd.m, dirpdf.m by Jan-Willem van de Meent (in tools/)
ML_loadStack2.m, by Fredrik Person (in +EMCCDfit/)
HMMcore/, by Martin Lindén
=====================================================================
---------------------------------------------------------------------
EMCCDfit contents
---------------------------------------------------------------------
localization_example.m  : example of spot localization
localization_example_framewise.m : spot localization frame by frame,
				   experimenting with different PSF
				   models and priors
MAP_EMCCD_refineSingleFrame.m 	 : localize spots in a single image
logL_EMCCD_lookup.m     : a lookup-table object for EMCCD noise
logL_psf.m              : a likelihood object that combines a camera noise 
                          model, a psf model, and an image
psf_diff_symgauss.m     : symmetric Gaussian PSF model
psf_diff_asymgauss_angle.m : asymmetric Gaussian PSF model, parameterized 
                             by 2 principal PSF widths and a rotation angle
log_likelihood_EMCCD_brute_parfor.m : Numerically compute EMCCD noise 
log_likelihood_EMCCD_brute_serial.m : likelihood, high gain approximation
ROItransform.m          : handle transformations in and out of subimages
ML_loadStack2.m         : read a tif stack
skewGauss_logPdf.m      : skew-normal probability density function
---------------------------------------------------------------------
EMhmm contents
---------------------------------------------------------------------
HMMexample.m			: example of a simple HMM analysis
EMhmm.preprocess.m		: preprocess SPT data
EMhmm.init_P_dat.m		: initialize an HMM
EMhmm.MLEmodelsearch.m		: high-level optimiztion for ML inference
EMhmm.MLEconverge.m		: converge a single HMM
EMhmm.positionEstimate.m 	: refined HMM-based estimate of localized
				  positions
EMhmm.diffusiveHMM_blur_detach.m: simulate SPT data from the
				  generative model of EMhmm

EMhmm.diffusionPathUpdate.m	: EM iteration to update diffusive path
EMhmm.hiddenStateUpdate.m	: EM iteration to update hidden states
EMhmm.MLEparameterUpdate.m	: EM iteration to update parameters
HMMcore/			: misc low-level HMM functions 
HMMcore/getDwellTRJ.m 		: extract dwell times etc from segmented 
				  hidden state trajectory
EMhmm.parameterEstimate.m	: Extract some parameter estimates
EMhmm.parameterBootstrap.m	: Bootstrap parameter estimates
EMhmm.positionEstimate.m	: Display parameter+-bootstrap err

----------------------------------------------------------------------
Installation: 
1) put the uncertainSPT folder and all subfiolders in the matlab path, e.g.
> cd path-to-uncertainSPT
> addpath(genpath(pwd))
Consult the Matlab documentation for ways to add the paths permanently.

2) recompile the HMMcore binaries if needed:
> cd path-to-uncertainSPT/HMMcore
> compile_code
---------------------------------------------------------------------
Release history

v0.9  (2016-07-08) : beta release for manuscript submission
v0.9.1(2016-08-03) : fixed repo error to include +EMCCDfit/ files
---------------------------------------------------------------------

