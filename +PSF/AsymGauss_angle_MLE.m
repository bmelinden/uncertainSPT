classdef AsymGauss_angle_MLE < PSF.AsymGauss_angle
    % An asymmetric Gaussian PSF model (see PSF.AsymGauss_angle) where the
    % background  intensity B has a log-normal prior (or log(B) is Normal)
    % with location parameter lnB0 and scale parameter lnBstd, which means
    % that <ln B> = lnB0, Std(ln B) lnBstd, <B> = exp(lnB0+lnBstd^2/2),
    % Std(B) = exp(lnB0+lnBstd^2/2)*sqrt(exp(lnBstd^2)-1),
    % see https://en.wikipedia.org/wiki/Log-normal_distribution.
    %
    % Construction:
    % P=PSF.AsymGauss_angle_MLE('initialGuess',[mux muy lnB lnN lnS1 lnS2 v])
    % P=PSF.AsymGauss_angle_MLE(...
    %     'initialGuess',[mux muy lnB lnN lnS1 lnS2 v],'priorParameters',[]);
    % (The initial guess is passed on to the AsymGauss_angle constructor.)
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = AsymGauss_angle_MLE(varargin)
            % call superclass constructor
            this@PSF.AsymGauss_angle(varargin{:});

            % sanity check on priors
            if(numel(this.priorParameters)>0)
               error('PSF.AsymGauss_angle_MLE needs empty or no priorParameters')
            end
        end
        
        % log prior distribution
        function [y,dy] = logPrior(this,param)
            y=0;
            dy=zeros(size(param));
        end
    end    
end