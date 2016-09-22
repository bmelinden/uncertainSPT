function param=psf2param_symgauss(p,param)
% param=psf2param_symgauss(p,param)
%
% Convert psf fit parameters to parameter struct, using a symmetric
% Gaussian PSF model, parameterized by 
% p=[ muX muY lnB lnN lnS] :
% E = exp(lnB) + f, with 
% f = N/2/pi/S^2*exp(-0.5*((xx-muX)/S)^2-0.5*((yy-muY)/S)^2)
%
% Note that this is a continuous model, so that exp(lnB) is the background
% intensity per unit area, not per pixel. But when doing math in pixel
% units the area per pixel is 1, and then this does not matter.
%
% param : input struct (optional). If given, existing fields are modofied
% (or not) accoridng to the output.
%
% output: parameter struct p, with fields
% p.background  : background intensity photons/area = exp(lnB)
% p.amplitude   : spot amplitude (photons)          = exp(lnN)
% p.std         : spot width, stdandard deviation   = exp(lnS)
%
% ML 2016-08-03

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMCCDfit.psf2param_symgauss.m, symmetric Gaussian PSF model parameters
% =========================================================================
% 
% Copyright (C) 2016 Martin Lindén
% 
% E-mail: bmelinden@gmail.com
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
if(~exist('param','var'))
    param=struct;
end
param.background=exp(p(3));
param.amplitude =exp(p(4));
param.std       =exp(p(5));
end
