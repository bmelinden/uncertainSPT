classdef SymGauss_logNormB < PSF.SymGauss
    % a symmetric Gaussian PSF model (SymGauss) where the background
    % intensity B has a log-normal prior with location parameter lnB0 and
    % scale parameter lnBstd, which means that
    % <ln B> = lnB0, Std(ln B) lnBstd, <B> = exp(lnB0+lnBstd^2/2),
    % Std(B) = exp(lnB0+lnBstd^2/2)*sqrt(exp(lnBstd^2)-1),
    % see https://en.wikipedia.org/wiki/Log-normal_distribution.
    %
    % Construction:
    % P = SymGauss_logNormB('initialGuess',[mux muy lnB lnN lnS],'priorParameters',[lnB0 lnBstd]),
    % (the initial guess is passed on to the SymGauss constructor)
    properties (Constant)
        priorName='log-normal background';
    end
    properties
        priorParameters	=[];
    end
    methods
        % Constructor
        function this = SymGauss_logNormB(varargin)
            % call superclass constructor
            this@PSF.SymGauss(varargin{:});
            % sanity check on priorParameters
            if(numel(this.priorParameters)~=2)
                error('PSF.SymGauss_logNormB needs priorParameters=[lnB0 lnBstd]')
            end
            if(this.priorParameters(2) <=0 )
                error('PSF.SymGauss_logNormB needs lnBstd >0')
            end
        end
        % compute log prior density
        function [y,dy] = logPrior(this,param)
           lnB=param(3);
           lnB0=this.priorParameters(1);
           lnBstd=this.priorParameters(2);
           
           y = - 1/2 *(lnB - lnB0).^2./(lnBstd^2);
           dy=zeros(size(param));
           dy(3)=- (lnB - lnB0)./(lnBstd^2);
        end        
    end
end