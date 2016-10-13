classdef SymGauss_logNormBN < PSF.SymGauss
    % a symmetric Gaussian PSF model (SymGauss) where the background
    % intensity B and spot amplitude Nhave log-normal priors: if B has
    % location parameter lnB0 and scale parameter lnBstd, 
    % then 
    % <ln B> = lnB0, Std(ln B) lnBstd, <B> = exp(lnB0+lnBstd^2/2),
    % Std(B) = exp(lnB0+lnBstd^2/2)*sqrt(exp(lnBstd^2)-1),
    % see https://en.wikipedia.org/wiki/Log-normal_distribution.
    %
    % Construction:
    % P = SymGauss_logNormB('initialGuess',[mux muy lnB lnN lnS],...
    %   'priorParameters',[lnB0 lnBstd lnN0 lnNstd]),
    % (the initial guess is passed on to the SymGauss constructor)
    properties (Constant)
        priorName='log-normal B,N';
    end
    properties
        priorParameters	=[];
    end
    methods
        % Constructor
        function this = SymGauss_logNormBN(varargin)
            % call superclass constructor
            this@PSF.SymGauss(varargin{:});
            % sanity check on priorParameters
            if(numel(this.priorParameters)~=4)
                error('PSF.SymGauss_logNormBN needs 4 priorParameters=[lnB0 lnBstd lnN0 lnNstd]')
            end
            if( ~isempty(find(this.priorParameters([2 4]) <=0,1)) )
                error('PSF.SymGauss_logNormBNS needs positive scale parameters lnBstd lnNstd.')
            end
        end
        % compute log prior density
        function [y,dy] = logPrior(this,param)
           lnB=param(3);
           lnB0=this.priorParameters(1);
           lnBvar=this.priorParameters(2)^2;
           
           lnN=param(4);
           lnN0=this.priorParameters(3);
           lnNvar=this.priorParameters(4)^2;

           %lnS=param(4);
           %lnS0=this.priorParameters(5);
           %Svar=this.priorParameters(6)^2;
           
           y = - 1/2 *(lnB - lnB0).^2./lnBvar -0.5*log(2*pi*lnBvar)...
               - 1/2 *(lnN - lnN0).^2./lnNvar -0.5*log(2*pi*lnNvar);%...
               %- 1/2 *(lnS - lnS0).^2./lnSvar -0.5*log(2*pi*lnSvar);
           dy=zeros(size(param));
           dy(3)=- (lnB - lnB0)./lnBvar;
           dy(4)=- (lnN - lnN0)./lnNvar;
           %dy(5)=- (lnS - lnS0)./lnSvar;
        end        
    end
end