function [E,dE_dmuX,dE_dmuY,dE_dlnB,dE_dlnN,dE_dlnS1,dE_dlnS2,dE_dv]=...
    psf_diff_asymgauss_angle(xx,yy,p)
% [E,dE_dmuX,dE_dmuY,dE_dlnN,dE_dlnS1,dE_dlnS2,dE_dv]=...
%           psf_diff_asymgauss_angle(xx,yy,p)
% asymmetric Gaussian PSF model and partial derivatives, 
% parameterized by parameters
% p=[ muX muY lnB lnN lnS1 lnS2 v] :
% E = exp(lnB) + f, with 
% f = N/2/pi/S1/S2*exp(-0.5*( cos(v)*(xx-muX)+sin(v)*(yy-muY))^2/S1^2
%                      -0.5*(-sin(v)*(xx-muX)+cos(v)*(yy-muY))^2/S2^2)

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


% ML: derivatives verified numerically 2916-05-31
muX=p(1);
muY=p(2);
lnB=p(3);
lnN=p(4);
lnS1=p(5);
lnS2=p(6);
v=p(7);

sig1=exp(lnS1);
sig2=exp(lnS2);
c =cos(v);
s =sin(v);
dx=xx-muX;
dy=yy-muY;

NEexp=1/2/pi*exp(lnN-lnS1-lnS2-1/2*(((c*dx+s*dy)/sig1).^2+((-s*dx+c*dy)/sig2).^2));
E=NEexp+exp(lnB);

if(nargout>1)
    dE_dmuX =NEexp.*( c/sig1^2*(c*dx+s*dy)-s/sig2^2*(-s*dx+c*dy));
    dE_dmuY =NEexp.*( s/sig1^2*(c*dx+s*dy)+c/sig2^2*(-s*dx+c*dy));
    dE_dlnB = ones(size(E))*exp(lnB);
    dE_dlnN = NEexp;
    dE_dlnS1 = NEexp.*(-1+((c*dx+s*dy)/sig1)^2);
    dE_dlnS2 = NEexp.*(-1+((-s*dx+c*dy)/sig2)^2);
    dE_dv = NEexp.*(c*dx+s*dy).*(s*dx-c*dy)*(1/sig1-1/sig2)*(1/sig1+1/sig2);
    % compute parameter dependent pixel intensities for symmetric Gaussian, with derivatives
end
