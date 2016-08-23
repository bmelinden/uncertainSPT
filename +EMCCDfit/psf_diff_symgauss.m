function [E,dE_dmuX,dE_dmuY,dE_dlnBG,dE_dlnN,dE_dlnS]=psf_diff_symgauss(xx,yy,p)
% [E,dE_dmuX,dE_dmuY,dE_dlnB,dE_dlnN,dE_dlnS]=psf_diff_psfGaussian(xx,yy,p)
%
% symmetric Gaussian PSF model and partial derivatives, parameterized by
% p=[ muX muY lnB lnN lnS] :
% E = exp(lnB) + f, with 
% f = N/2/pi/S^2*exp(-0.5*((xx-muX)/S)^2-0.5*((yy-muY)/S)^2)
%
% Note that this is a continuous model, so that exp(lnB) is the background
% intensity per unit area, not per pixel. But when doing math in pixel
% units the area per pixel is 1, and then this does not matter.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMCCDfit.psf_diff_symgauss.m, symmetric Gaussian PSF model
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


muX=p(1);
muY=p(2);
lnB=p(3);
N=exp(p(4));
S2=exp(2*p(5)); % S^2

% 2D Gaussian density profile
NEexp=N/S2*exp(-1/2/S2*((muX-xx).^2+(muY-yy).^2))/2/pi;

E=exp(lnB)+NEexp;

if(nargout>1)
    dE_dmuX =-(muX-xx).*NEexp/S2;
    dE_dmuY =-(muY-yy).*NEexp/S2;
    dE_dlnBG=ones(size(E))*exp(lnB);
    dE_dlnN = NEexp;
    dE_dlnS = NEexp.*(-2+((muX-xx).^2+(muY-yy).^2)/S2);
    % compute parameter dependent pixel intensities for symmetric Gaussian, with derivatives
    % ML 2015-11-17 : partial derivatives validated
end
end
