classdef SymGaussS0_logNormN_expS0 < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0), with
    % an exponential prior on the rescaled psf width y=(S-S0)/S0=exp(lnS),
    % which translates to a prior density  
    %
    % p0(lnS) = 1/y0*exp(lnS-exp(lnS)/y0)
    % 
    % for an exponential prior with mean value y0. The spot amplitude has a
    % log-normal priors, i.e., a normal prior on the fit parameter lnN,
    % with mean value parameter lnN0 and std-parameter lnNstd.
    %
    % Construction:
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lnS],'S0','priorParameters',[lnN0 lnNstd y0]),
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lnS],'lambda',lambda,'NA',NA,'priorParameters',[lnB0 lnBstd y0]),
    % (the initial guess, NA,lambda or S0 are passed on to the SymGaussS0 constructor)
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = SymGaussS0_logNormN_expS0(varargin)
            % call superclass constructor
            this@PSF.SymGaussS0(varargin{:});
            % sanity check on priors
            if(numel(this.priorParameters)~=3 )
               error('PSF.SymGaussS0_logNormN_expS0 needs 3 prior parameters [lnN0 lnNstd y0]')
            end
        end
        % log prior density
        function [y,dy] = logPrior(~,param)
            
            lnN=param(4);
            lnN0=this.priorParameters(1);
            lnNvar=this.priorParameters(2)^2;
           
            y0=this.priorParameters(3);            
            
            y= -1/2 *(lnN - lnN0).^2./lnNvar -0.5*log(2*pi*lnNvar)...
               -log(y0)+lnS-exp(lnS)/y0;
           
            dy=zeros(size(param));
            dy(4)=- (lnN - lnN0)./lnNvar;
            dy(5) = 1-exp(lnS)/y0;
        end
    end
end