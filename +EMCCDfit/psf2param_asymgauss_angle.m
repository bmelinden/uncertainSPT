function param=psf2param_asymgauss_angle(p,param)
% param=psf2param_asymgauss_angle(p,param)
% asymmetric Gaussian PSF model and partial derivatives, 
% parameterized by parameters
% p=[ muX muY lnB lnN lnS1 lnS2 v] :
% E = exp(lnB) + f, with 
% f = N/2/pi/S1/S2*exp(-0.5*( cos(v)*(xx-muX)+sin(v)*(yy-muY))^2/S1^2
%                      -0.5*(-sin(v)*(xx-muX)+cos(v)*(yy-muY))^2/S2^2)
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
% p.std1        : spot principal width 1, standard deviation   = exp(lnS1)
% p.std2        : spot principal width 2, standard deviation   = exp(lnS2)
% p.angle       : spot orientation, = v
% ML 2016-08-03

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMCCDfit.psf_diff_asymgauss.m, asymmetric Gaussian PSF model
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
lnB=p(3);
lnN=p(4);
lnS1=p(5);
lnS2=p(6);
v=p(7);

sig1=exp(lnS1);
sig2=exp(lnS2);

if(~exist('param','var'))
    param=struct;
end
param.background=exp(lnB);
param.amplitude =exp(lnN);
param.std1      =exp(lnS1);
param.std2      =exp(lnS2);
param.angle     =v;

