classdef SymGaussS0_logNormBN_expS0 < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0), with
    % an exponential prior on the excess psf width dS=(S-S0)/S0=exp(lndS),
    % which translates to a prior density  
    %
    % p0(lndS) = 1/dS0*exp(lndS-exp(lndS)/dS0)
    % 
    % for an exponential prior with mean value dS0. 
    %
    % The background and spot amplitude have log-normal priors, i.e.,
    % normal priors on the fit parameters lnB, lnN, with mean value
    % parameters lnB0, lnN0, and std-parameters lnBstd, lnNstd. Setting
    % lnBstd, lnNstd to inf gives flat priors on lnB or lnN, respectively.
    % lnBstd=lnNstd=inf is equivalent to SymGaussS0_expS0.m
    %
    % Construction:
    % P = SymGaussS0_logNormBN_expS0('initialGuess',[mux muy lnB lnN lndS],'S0','priorParameters',[lnB0 lnBstd lnN0 lnNstd dS0]),
    % P = SymGaussS0_logNormBN_expS0('initialGuess',[mux muy lnB lnN lndS],'lambda',lambda,'NA',NA,'priorParameters',[lnB0 lnBstd lnN0 lnNstd dS0]),
    % (the initial guess, NA,lambda or S0 are passed on to the SymGaussS0 constructor)
    properties (Constant)
        priorName='logNorm(B,N), exp(dS)';
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
               error('PSF.SymGaussS0_logNormN_expS0 needs 5 prior parameters [lnB0 lnBstd lnN0 lnNstd dS0]')
            end
        end
        % log prior density
        function [y,dy] = logPrior(this,param)
            
            lndS=param(5);
            dS0=this.priorParameters(5);
            
            y= -log(dS0)+lndS-exp(lndS)/dS0;
            dy=zeros(size(param));
            dy(5) = 1-exp(lndS)/dS0;
            
            lnBvar=this.priorParameters(2)^2;
            if(isfinite(lnBvar))
                lnB=param(3);
                lnB0=this.priorParameters(1);
                y=y -1/2 *(lnB - lnB0).^2./lnBvar -0.5*log(2*pi*lnBvar);
                dy(3)=- (lnB - lnB0)./lnBvar;
            end
            
            lnNvar=this.priorParameters(4)^2;
            if(isfinite(lnNvar))
                lnN=param(4);
                lnN0=this.priorParameters(3);
                y=y -1/2 *(lnN - lnN0).^2./lnNvar -0.5*log(2*pi*lnNvar);
                dy(4)=- (lnN - lnN0)./lnNvar;
            end
        end
    end
end
