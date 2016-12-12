Diffusion constant estimator for 1-state diffusive SPT data with 
localization errors and blur.

Dest.logLlambda - log likelihood function for camera-based single
particle tracking data with camera-based blur and localization errors
that are given for every point, similar to those described in [1].

This code makes use of the routine triSym_d1Inv_trjWise.mex (source
code in HMMcore) to compute the determinant and some elements of the
inverse of a symmetric tridiagonal matrix. This is a version of the Thomas 
algorithm (fast, but may become unstable in some situations).

Dest.preprocess_mixed_columns : data preprocessor.

Note: if you want to look under the hood, note that these functions may use 
partly different index conventions than those in +EMhmm, and expect 
precision as STANDARD DEVIATIONS, while EMhmm work with variances.

Martin Lindén, bmelinden@gmail.com, 2016-12-12

1. P. K. Relich, M. J. Olah, P. J. Cutler, and K. A. Lidke, "Estimation of 
the diffusion constant from intermittent trajectories with variable 
position uncertainties," Phys. Rev. E, vol. 93, no. 4, p.042401, Apr. 2016.
