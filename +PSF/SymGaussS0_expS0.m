classdef SymGaussS0_expS0 < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0), with
    % an exponential prior on the excess psf width 
    % dS = (S-S0)/S0 = exp(lndS), which translates to a prior density 
    %
    % p0(lndS) = 1/dS0*exp(lndS-exp(lndS)/dS0), for an exponential prior with mean
    % value <dS>=dS0
    %
    % Construction:
    % P = SymGaussS0_expS0('initialGuess',[mux muy lnB lnN lndS],'S0','priorParameters',dS0),
    % P = SymGaussS0_expS0('initialGuess',[mux muy lnB lnN lndS],'lambda',lambda,'NA',NA,'priorParameters',dS0),
    % (the initial guess, NA,lambda,S0 are passed on to the SymGaussS0 constructor)
    properties (Constant)
        priorName='exp(dS)';
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
               error('PSF.SymGaussS0_exp needs a single non-negative prior parameter dS0')
            end
        end
        % log prior density
        function [y,dy] = logPrior(this,param)
            
            lndS=param(5);
            dS0=this.priorParameters(1);            
            
            y=-log(dS0)+lndS-exp(lndS)/dS0;
            dy=zeros(size(param));
            dy(5) = 1-exp(lndS)/dS0;
        end
    end
end
