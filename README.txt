---------------------------------------------------------------------
uncertainSPT
---------------------------------------------------------------------
Uncertain SPT contains two software tools to extract and use
localization uncertainty for localization microscopy and single
particle tracking.
1) EMCCDfit, a localization algorithm to estimate particle positions
and from images aquired by an EMCCD camera, and
2) EMhmm, a variational EM algorithm that performs maximum likelihood
(ML) inference in a diffusive hidden Markov model, which accounts for
motion blur and localization uncertainty.
======================================================================
Copyright (C) 2016 Martin Lindén and Johan Elf

E-mail: bmelinden@gmail.com, johan.elf@gmail.com
======================================================================
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or any
later version.  This program is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

If you use vbSPT please cite: "Person F, Lindén M, Unoson C, Elf J,
Extracting intracellular reaction rates from single molecule tracking
data, Nature Methods 10, 265–269 (2013). doi:10.1038/nmeth.2367"

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

This product includes code developed by Jan-Willem van de Meent (in
external/), and Fredrik Persson (in +EMCCDfit/), see copyright in
individual files.
=====================================================================
---------------------------------------------------------------------
EMCCDfit contents
---------------------------------------------------------------------
localization_example.m  : example of spot localization
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
----------------------------------------------------------------------
Installation: put the uncertainSPT folder and all subfiolders in the
matlab path, e.g.
> cd path-to-uncertainSPT
> addpath(genpath(pwd))
(consult the Matlab documentation for ways to make it permanent).
---------------------------------------------------------------------
Release history

v1.0  (2016-07-xx) : initial release


---------------------------------------------------------------------
to do list and notes (to be removed before first release)

- fix names etc in  
+EMhmm/preprocess.m
+EMhmm/diffusiveHMM_blur_detach.m

- add license information to all files
- create a simple HMM example
- create a brief spot-fit example with a short movie
- create a brief documentation: functions, and what they do
- do I need log_likelihood_EMCCD_brute_parfor ? Or EMCCDfit.logL_EMCCD_lookup?

copyright suggestion:
%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB3_switching_diffusion_ddim.m, simulate d-dimensional diffusion with 
% multiple diffusion constants, part of the vbSPT package
% =========================================================================
% 
% Copyright (C) 2013 Martin Lindén and Fredrik Persson
% 
% E-mail: bmelinden@gmail.com, freddie.persson@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by linking or combining it
%  with Matlab or any Matlab toolbox, the licensors of this Program grant you 
%  additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% start of actual code
