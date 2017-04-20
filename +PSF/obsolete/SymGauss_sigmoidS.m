classdef SymGauss_sigmoidS < PSF.SymGauss
    % a symmetric Gaussian PSF model (SymGauss) with an improper sigmoid
    % prior on the PSF width,
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
    % P = SymGauss_MLE('initialGuess',[mux muy lnB lnN lnS],'priorParameters',[dsScale S0]),
    % P = SymGauss_MLE('initialGuess',[mux muy lnB lnN lnS],'priorParameters',[dsScale lambda NA]),
    % (the initial guess is passed on to the SymGauss constructor)
    properties (Constant)
        priorName='sigmoid(S,sScale)';
    end
    properties
        priorParameters	=[];
        lambda=[];
        NA=[];
    end
    
    methods
        % Constructor
        function this = SymGauss_sigmoidS(varargin)
            % call superclass constructor
            this@PSF.SymGauss(varargin{:});
            % interpret priors
            if(numel(this.priorParameters)==2)     % then [dsScale S0]: do nothing
            elseif(numel(this.priorParameters)==3) % then [dsScale lambda NA]: compute S0
                this.lambda=this.priorParameters(2);
                this.NA    =this.priorParameters(3);
                
                S0=0.21*this.lambda/this.NA;
                this.priorParameters=[this.priorParameters(1) S0];
            else
               error('PSF.SymGauss_sigmoidS needs 2 priorParameters [S0 dsScale]')
            end
        end
        % log prior density
        function [y,dy] = logPrior(this,param)
            
            dsScale=this.priorParameters(1);
            S0=this.priorParameters(2);
            lnS=param(5);
            S=exp(lnS);
            
            y=-log1p(exp(-(S/S0-1)/dsScale));
            dy=zeros(size(param));
            dy(5)=S/dsScale/S0./(1+exp((S/S0-1)/dsScale));
        end
    end
end