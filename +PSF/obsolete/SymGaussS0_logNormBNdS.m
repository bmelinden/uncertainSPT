classdef SymGaussS0_logNormBNdS < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0), with
    % an exponential prior on the excess psf width dS=(S-S0)/S0=exp(lndS)
    %
    % Priors: log-normal on all parameters, i.e., normal priors on the
    % fit parameters lnB, lnN, lndS, with mean value parameters lnN0, lnB0,
    % lndS0, and std-parameters lnNstd, lnBstd, lndSstd. Setting a
    % std-parameter to inf leads to a flat prior on the corresponding
    % parameter. lnNstd=lnBstd=lndSstd=inf is equivalent to SymGaussS0_MLE.
    %
    % Construction:
    % P = SymGaussS0_logNormBNdS('initialGuess',[mux muy lnB lnN lndS],'S0','priorParameters',[lnB0 lnBstd lnN0 lnNstd lndS0 lndSstd]), 
    % P = SymGaussS0_logNormBNdS('initialGuess',[mux muy lnB lnN lndS],'lambda',lambda,'NA',NA,'priorParameters',[lnB0 lnBstd lnN0 lnNstdlndS0 lndSstd]),
    % (the initial guess, NA,lambda or S0 are passed on to the SymGaussS0
    % constructor).
    properties (Constant)
        priorName='logNorm(B,N,dS)';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = SymGaussS0_logNormBNdS(varargin)
            % call superclass constructor
            this@PSF.SymGaussS0(varargin{:});
            % sanity check on priors
            if(numel(this.priorParameters)~=4 )
               error('PSF.SymGaussS0_logNormBNdS needs 4 prior parameters [lnB0 lnBstd lnN0 lnNstd]')
            end
        end
        % log prior density
        function [y,dy] = logPrior(this,param)
            
            y=0;
            dy=zeros(size(param));
            
            lnBvar=this.priorParameters(2)^2;
            lnNvar=this.priorParameters(4)^2;
            lndSvar=this.priorParameters(6)^2;

            if(isfinite(lnBvar))
                lnB=param(3);
                lnB0=this.priorParameters(1);
                y=y -1/2 *(lnB - lnB0).^2./lnBvar -0.5*log(2*pi*lnBvar);
                dy(3)=- (lnB - lnB0)./lnBvar;
            end
            
            if(isfinite(lnNvar))
                lnN=param(4);
                lnN0=this.priorParameters(3);
                y=y -1/2 *(lnN - lnN0).^2./lnNvar -0.5*log(2*pi*lnNvar);
                dy(4)=- (lnN - lnN0)./lnNvar;
            end

            if(isfinite(lndSvar))
                lndS=param(5);
                lnB0=this.priorParameters(5);
                y=y -1/2 *(lndS - lndS0).^2./lndSvar -0.5*log(2*pi*lndSvar);
                dy(5)=- (lndS - lndS0)./lndSvar;
            end
        end
    end
end
