classdef SymGaussS0_logNormBN_expS0 < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0), with
    % an exponential prior on the rescaled psf width y=(S-S0)/S0=exp(lnS),
    % which translates to a prior density  
    %
    % p0(lnS) = 1/y0*exp(lnS-exp(lnS)/y0)
    % 
    % for an exponential prior with mean value y0. The background and spot
    % amplitude have log-normal priors, i.e., normal priors on the fit
    % parameters lnB, lnN, with mean value parameters lnB0, lnN0, and
    % std-parameters lnBstd, lnNstd.
    %
    % Construction:
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lnS],'S0','priorParameters',[lnB0 lnBstd lnN0 lnNstd y0]),
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lnS],'lambda',lambda,'NA',NA,'priorParameters',[lnB0 lnBstd lnN0 lnNstd y0]),
    % (the initial guess, NA,lambda or S0 are passed on to the SymGaussS0 constructor)
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = SymGaussS0_logNormBN_expS0(varargin)
            % call superclass constructor
            this@PSF.SymGaussS0(varargin{:});
            % sanity check on priors
            if(numel(this.priorParameters)~=5 )
               error('PSF.SymGaussS0_logNormN_expS0 needs 5 prior parameters [lnB0 lnBstd lnN0 lnNstd y0]')
            end
        end
        % log prior density
        function [y,dy] = logPrior(~,param)
            
            lnB=param(3);
            lnB0=this.priorParameters(1);
            lnBvar=this.priorParameters(2)^2;
           
            lnN=param(4);
            lnN0=this.priorParameters(3);
            lnNvar=this.priorParameters(4)^2;
           
            y0=this.priorParameters(5);            
            
            y= -1/2 *(lnB - lnB0).^2./lnBvar -0.5*log(2*pi*lnBvar)...
               -1/2 *(lnN - lnN0).^2./lnNvar -0.5*log(2*pi*lnNvar)...
               -log(y0)+lnS-exp(lnS)/y0;
           
           dy=zeros(size(param));
           dy(3)=- (lnB - lnB0)./lnBvar;
           dy(4)=- (lnN - lnN0)./lnNvar;
           dy(5) = 1-exp(lnS)/y0;
        end
    end
end