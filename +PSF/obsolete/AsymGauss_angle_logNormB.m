classdef AsymGauss_angle_logNormB < PSF.AsymGauss_angle
    % An asymmetric Gaussian PSF model (see AsymGauss_angle) where the
    % background  intensity B has a log-normal prior (or log(B) is Normal)
    % with location parameter lnB0 and scale parameter lnBstd, which means
    % that <ln B> = lnB0, Std(ln B) lnBstd, <B> = exp(lnB0+lnBstd^2/2),
    % Std(B) = exp(lnB0+lnBstd^2/2)*sqrt(exp(lnBstd^2)-1),
    % see https://en.wikipedia.org/wiki/Log-normal_distribution.
    %
    % Construction:
    % P = SymGauss_logNormB(...
    %   'priorParameters',[lnB0 lnBstd],'initialGuess',[mux muy lnB lnN lnS]
    % (The initial guess is passed on to the AsymGauss constructor.)
    properties (Constant)
        priorName='log-normal background';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = AsymGauss_angle_logNormB(varargin)
            % call superclass constructor
            this@PSF.AsymGauss_angle(varargin{:});
            % sanity check on prior parameters
            if(numel(this.priorParameters)~=2)
                error('PSF.AsymGauss_angle needs priorParameters=[lnB0 lnBstd]')
            end
            if(this.priorParameters(2) <=0 )
                error('PSF.AsymGauss_angle needs lnBstd >0')
            end
        end
        
        % log prior distribution
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
