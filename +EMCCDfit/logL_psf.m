%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMCCDfit.logL_psf, log likelihood lookup table for high gain EMCCD noise
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

classdef logL_psf
    properties
        logLobj=[];
        psfFun=[];
        
        % the image
        image=[];
        imX=[];
        imY=[];
        imSize=[];
        
        % quadrature quantities
        Xquad=[];
        Yquad=[];
        W=[];
    end
    methods
        % constructor
        function obj=logL_psf(logL,psfFun,C,Nquad)
            % Create a fast dot fit object.
            % logL  : a log-likelihood object, with method
            %         logL.lnL(image,intensity), e.g., an instance of
            %         EMCCDfit.logL_EMCCD_lookup.m
            % psfFun: a psf model function handle, of the form
            %         E=psfFun(x,y,p),
            %         where p are all psf model parameters
            % C     : an image, consistent with the logL camera model
            % Nquad : number of quadrature points per pixel in each
            %         dimension (i.e., there will be Nquad^2 psfFun
            %         evaluations per pixel). Quadrature points are equally
            %         spaced.
            %
            % Methods
            %
            % lnL=logL_psf(obj,p) : log likelihood for parameters p (as
            % defined by psfFun).  NOTE: many standard optimization
            % routines such as fminunc does minimization, which meanns the
            % objective function needs to be -lnL.
            %
            % [E,C]=psfModel(obj,p) : evaluate the PSF intensity image E
            % corresponding to parameters p (as above), and return the
            % image data C.
            % ML 2016-05-31
            
            % store a reference to the log likelohood function
            obj.logLobj=logL;
            
            % store the psf function handle
            obj.psfFun=psfFun;
            
            % store image data
            obj.imSize=size(C);
            rows=1:size(C,1);
            cols=1:size(C,2);
            [X,Y]=meshgrid(cols,rows);
            obj.image=reshape(C,numel(C),1);
            obj.imX  =reshape(X,numel(X),1);
            obj.imY  =reshape(Y,numel(Y),1);
            
            % precompute pixel integration weights and
            qPts=linspace(-0.5,0.5,2*Nquad+1);
            qPts=qPts(2:2:end-1);
            [dx,dy]=meshgrid(qPts,qPts);
            dx=reshape(dx,1,Nquad^2);
            dy=reshape(dy,1,Nquad^2);
            
            % grid for pixel quadrature
            [a,b]=meshgrid(dx,obj.imX);
            obj.Xquad=a+b;
            [a,b]=meshgrid(dy,obj.imY);
            obj.Yquad=a+b;
            obj.W=1/Nquad^2; % same weight for all quadrature points here
        end
        % parameter model
        function [E,C]=psfModel(obj,p)
            % compute the model corresponding to parameters p
            %
            % E: PSF model (units photons/pixel)
            % C: image data (counts)

            E=obj.psfFun(obj.Xquad,obj.Yquad,p);
            E=reshape(sum(E,2)*obj.W,obj.imSize(1),obj.imSize(2));
            C=reshape(obj.image,obj.imSize(1),obj.imSize(2));
        end
        % model likelihood 
        function L=lnL(obj,p)
            % compute log-likelihood corresponding to
            % parameters p.

            E=obj.psfFun(obj.Xquad,obj.Yquad,p);
            % quadrature over single pixels
            E=sum(E,2)*obj.W;
            % compute likelihood
            lnLi=obj.logLobj.lnL(obj.image,E);
            L=sum(lnLi);
        end
    end
end
