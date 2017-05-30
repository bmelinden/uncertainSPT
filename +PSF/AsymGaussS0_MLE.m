classdef AsymGaussS0_MLE < PSF.AsymGaussS0
    % An asymmetric Gaussian PSF model with PSF width offset (see
    % PSF.AsymGaussS0) and flat prior on the fit parameters
    % lnB,lnN,lndS1,lndS2,v
    %
    % Construction:
    % P=PSF.AsymGauss_angle_MLE('initialGuess',[mux muy lnB lnN lndS1 lndS2 v])
    % P=PSF.AsymGauss_angle_MLE(...
    %     'initialGuess',[mux muy lnB lnN lndS1 lndS2 v],'priorParameters',[]);
    % (The initial guess is passed on to the AsymGauss_lndS_angle constructor.)
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = AsymGaussS0_MLE(varargin)
            % call superclass constructor
            this@PSF.AsymGaussS0(varargin{:});

            % sanity check on priors
            if(numel(this.priorParameters)>0)
               error('PSF.AsymGaussS0_MLE needs empty or no priorParameters')
            end
        end
        
        % log prior distribution
        function [y,dy] = logPrior(this,param)
            y=0;
            dy=zeros(size(param));
        end
    end    
end