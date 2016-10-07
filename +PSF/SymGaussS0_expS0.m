classdef SymGaussS0_expS0 < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0), with
    % an exponential prior on the rescaled psf width 
    % y = (S-S0)/S0 = exp(lnS), which translates to a prior density 
    %
    % p0(lnS) = 1/y0*exp(lnS-exp(lnS)/y0), for an exponential prior with 
    %
    % Construction:
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lnS],'S0','priorParameters',y0),
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lnS],'lambda',lambda,'NA',NA,'priorParameters',y0),
    % (the initial guess, NA,lambda,S0 are passed on to the SymGaussS0 constructor)
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = SymGaussS0_expS0(varargin)
            % call superclass constructor
            this@PSF.SymGaussS0(varargin{:});
            % sanity check on priors
            if(numel(this.priorParameters)~=1 || this.priorParameters <=0)
               error('PSF.SymGaussS0_exp needs a single non-negative prior parameter')
            end
        end
        % log prior density
        function [y,dy] = logPrior(~,param)
            y0=this.priorParameters(1);            
            y=-log(y0)+lnS-exp(lnS)/y0;
            dy=zeros(size(param));
            dy(5) = 1-exp(lnS)/y0;
        end
    end
end