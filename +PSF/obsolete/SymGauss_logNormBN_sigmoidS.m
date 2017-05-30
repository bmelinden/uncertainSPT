classdef SymGauss_logNormBN_sigmoidS < PSF.SymGauss
    % a symmetric Gaussian PSF model (SymGauss) where the background
    % intensity B and spot amplitude N have log-normal priors with location
    % parameters lnB0, lnN0, and scale parameters lnBstd, lnNstd. Choosing
    % a scale parameter = inf turns off the corresponding prior.
    % 
    % The PSF width takes an improper sigmoid prior, 
    % p(S)  = 1/(1+exp(-(S/S0-1)/dsScale))/S, or
    % p(lnS)= 1/(1+exp(-(exp(lnS)/S0-1)/dsScale))
    %
    % This effectively enforces a lower bound on S, which can be chosen to
    % be that of a well-focused spot, S0 = 0.21*lambda/NA, and so the prior
    % can be constructed in two possible ways.
    % The factor 1/S is needed in order that the prior on log(S) must not
    % grow as S>>S0, which could potentially cause lnS to grow without
    % bound during the optimization. 
    %
    % Construction:
    % P = SymGauss_logNormB_sigmoidS('initialGuess',[mux muy lnB lnN lnS],'priorParameters',[lnB0 lnBstd lnN0 lnNstd dsScale S0]),
    % P = SymGauss_logNormB_sigmoidS('initialGuess',[mux muy lnB lnN lnS],'priorParameters',[lnB0 lnBstd lnN0 lnNstd dsScale lambda NA]),
    % (the initial guess is passed on to the SymGauss constructor)
    properties (Constant)
        priorName='log-normal B,N, sigmoid S';
    end
    properties
        priorParameters	=[];
        lambda=[];
        NA=[];
    end
    methods
        % Constructor
        function this = SymGauss_logNormBN_sigmoidS(varargin)
            % call superclass constructor
            this@PSF.SymGauss(varargin{:});
            % sanity check on priorParameters
            if(numel(this.priorParameters)==6 ) % then S0 is specified
            elseif(numel(this.priorParameters)==7 ) % then compute S0
                this.lambda=this.priorParameters(6);
                this.NA    =this.priorParameters(7);
                
                S0=0.21*this.lambda/this.NA;
                this.priorParameters=[this.priorParameters(1:5) S0];
            else
                error('PSF.SymGauss_logNormBNS_sigmoidS needs priorParameters=[lnB0 lnBstd lnN0 lnNstd dsScale S0] or [lnB0 lnBstd lnN0 lnNstd dsScale lambda NA]')
            end
        end
        % compute log prior density
        function [y,dy] = logPrior(this,param)
            
            dsScale=this.priorParameters(5);
            S0=this.priorParameters(6);
            lnS=param(5);
            S=exp(lnS);
            
            y=-log1p(exp(-(S/S0-1)/dsScale));
            dy=zeros(size(param));
            dy(5)=S/dsScale/S0./(1+exp((S/S0-1)/dsScale));
            
            lnBvar=this.priorParameters(2)^2;
            if(isfinite(lnBvar))
                lnB=param(3);
                lnB0=this.priorParameters(1);
                y = y - 1/2 *(lnB - lnB0).^2./lnBvar -0.5*log(2*pi*lnBvar);
                dy(3)=- (lnB - lnB0)./lnBvar;
            end
            
            lnNvar=this.priorParameters(4)^2;
            if(isfinite(lnNvar))
                lnN=param(4);
                lnN0=this.priorParameters(3);
                y=y - 1/2 *(lnN - lnN0).^2./lnNvar -0.5*log(2*pi*lnNvar);
                dy(4)=- (lnN - lnN0)./lnNvar;
            end
        end
    end
end
