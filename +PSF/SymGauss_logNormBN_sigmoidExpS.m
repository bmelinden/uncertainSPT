classdef SymGauss_logNormBN_sigmoidExpS < PSF.SymGauss
    % a symmetric Gaussian PSF model (SymGauss) with 
    % - log-normal prior distributions on background and spot amplitudes,
    % with mean value- and scale parameters lnB0, lnN0, and lnBstd, lnNstd,
    % respectively. lnBstd, lnNstd=inf means that no prior is applied.
    % - a softly truncated exponential prior on the relative PSF width
    % s=S/S0, where S0 is a lower bound on the spot width S. The truncated
    % exponential prior would be that s-1|s>1 ~ exp(L), i.e., 
    %
    % p0(s) = 1/L*exp(-s/L)/(1-exp(-1/L)), s>=1, 
    %
    % but instead of enforcing this lower bound explicitly, we use a soft
    % truncation with a sigmoid function 
    %
    % b(s) = 1/(1+ exp(-(s-1)/dsScale)).
    %
    % If dsScale << L, the distribution b(s)*p0(s) is approximately
    % normalized (and for parameter MAP optimization, the normalization
    % constant is not relevant anyway), and the prior on the actual fit
    % parameter lnS=log(s)+log(S0) becomes 
    %
    % p0(lnS) = exp(-1/L*(exp(lnS)/S0-1)) / 
    %           (1 + exp(-1/dsScale*(exp(lnS)/S0-1))).
    %
    % This effectively enforces a S >= S0, and a physically reasonable
    % choise is the width of a well-focused spot, S0 = 0.21*lambda/NA. We
    % therefore allow two different prior constructions:
    %
    % P = SymGauss_MLE('initialGuess',[mux muy lnB lnN lnS],'priorParameters',[lnB0 lnBstd lnN0 lnNstd L dsScale S0]),
    % P = SymGauss_MLE('initialGuess',[mux muy lnB lnN lnS],'priorParameters',[lnB0 lnBstd lnN0 lnNstd L dsScale lambda NA]),
    %
    % ML 2016-10-13
    
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
        lambda=[];
        NA=[];
    end
    
    methods
        % Constructor
        function this = SymGauss_logNormBN_sigmoidExpS(varargin)
            % call superclass constructor
            this@PSF.SymGauss(varargin{:});
            % interpret priors
            if(numel(this.priorParameters)==7)     % then [lnB0 lnBstd lnN0 lnNstd L dsScale S0]: do nothing
            elseif(numel(this.priorParameters)==8) % then [lnB0 lnBstd lnN0 lnNstd L dsScale lambda NA]: compute S0
                
                this.lambda=this.priorParameters(6);
                this.NA    =this.priorParameters(7);
                
                S0=0.21*this.lambda/this.NA;
                this.priorParameters=[this.priorParameters(1:6) S0];
            else
               error('PSF.SymGauss_sigmoidS needs priorParameters [L S0 dsScale] or [L dsScale lambda NA]')
            end
        end
        % log prior density
        function [y,dy] = logPrior(this,param)
            
            lnB    =param(3);
            lnB0   =this.priorParameters(1);
            lnBvar =this.priorParameters(2)^2;
            
            lnN    =param(4);
            lnN0   =this.priorParameters(3);
            lnNvar =this.priorParameters(4)^2;
            
            lnS=param(5);
            S=exp(lnS);
            L      =this.priorParameters(5);
            dsScale=this.priorParameters(6);
            S0     =this.priorParameters(7);
            
            
            y = -log(L)-(S/S0-1)/L-log1p(exp(-(S/S0-1)/dsScale));
            dy=zeros(size(param));
            dy(5)=S/dsScale/S0./(1+exp((S/S0-1)/dsScale))-S/L/S0;
            % derivative numerically verified, ML 2016-10-10
            
            if(isfinite(lnBvar))
                y=y-1/2 *(lnB - lnB0).^2./lnBvar -0.5*log(2*pi*lnBvar);
                dy(3)=- (lnB - lnB0)./lnBvar;
            end
            if(isfinite(lnNvar))
                y=y-1/2 *(lnN - lnN0).^2./lnNvar -0.5*log(2*pi*lnNvar);
                dy(4)=- (lnN - lnN0)./lnNvar;
            end
            
        end
    end
end