function [rowInd,colInd,x0,y0]=ROItransform(x_c,y_c,x_w,y_w,imSize)
% [rowInd,colInd,x0,y0]=ROItransform(x_c,y_c,x_w,y_w,imSize)
%
% Compute index range and coordinate transform for a region of interest
% (ROI) in an image.
%
% input
% x_c,y_c   : ROI center, in pixel coordinates
% x_w,y_w   : ROI width in the x and y directions (pixels)
% imSize    : Size of the parent image, as in imSize=size(parentImage).
%             (Only the first two indices are used).
%
% Output
% rowInd,colInd   : The region of interest can be plotted by, e.g., 
%                   imagesc(parentImage(rowInd,colInd). The index ranges
%                   are truncated to be in the allowed range of the parent
%                   image.
% x0,y0 : coordinate transform between parent and ROI.
%         x_parent = x_ROI+x0, y_parent=y_ROI+y0.
%
% Martin Lindén, bmelinden@gmail.com, 2015-11-05

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMCCDfit.ROItransform.m, indices to sub-images for localization
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


x_min=max(1        , ceil(x_c-x_w/2));
x_max=min(imSize(2),floor(x_c+x_w/2));

y_min=max(1        , ceil(y_c-y_w/2));
y_max=min(imSize(1),floor(y_c+y_w/2));

rowInd=y_min:y_max;
colInd=x_min:x_max;

x0=x_min-1;
y0=y_min-1;


