classdef AsymGauss_angle_logNormBNS < PSF.AsymGauss_angle
    % An asymmetric Gaussian PSF model (see AsymGauss_angle) where all
    % shape parameters (background B, amplitude N, psf widths S1,S2) except
    % the asymmetry orientation angle have log-normal priors. A log-normal
    % distribution with location parameter lnB0 and scale parameter lnBstd,
    % has 
    % <ln B> = lnB0, Std(ln B) lnBstd, <B> = exp(lnB0+lnBstd^2/2),
    % Std(B) = exp(lnB0+lnBstd^2/2)*sqrt(exp(lnBstd^2)-1),
    % see https://en.wikipedia.org/wiki/Log-normal_distribution.
    %
    % Construction:
    % P = SymGauss_logNormBNS(...
    %   'priorParameters',[lnB0 lnBstd lnN0 lnNstd lnS0 lnSstd],...
    %   'initialGuess',[mux muy lnB lnN lnS1 lnS2 v]
    % (The initial guess is passed on to the AsymGauss constructor.)
    % Note that both psf width S1,S2 have the same prior parameters
    properties (Constant)
        priorName='log-normal background';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = AsymGauss_angle_logNormBNS(varargin)
            % call superclass constructor
            this@PSF.AsymGauss_angle(varargin{:});
            % sanity check on prior parameters
            if(numel(this.priorParameters)~=6)
                error('PSF.AsymGauss_angle_logNormBNS needs priorParameters=[lnB0 lnBstd lnN0 lnNstd lnS0 lnSstd]')
            end
            if( ~isempty(find(this.priorParameters(2:2:6) <=0,1)) )
                error('PSF.AsymGauss_angle_logNormBNS needs positive scale parameters lnBstd lnNstd lnSstd.')
            end
        end
        
        % log prior distribution
        function [y,dy] = logPrior(this,param)

           lnB=param(3);
           lnN=param(4);
           lnS1=param(5);
           lnS2=param(6);

           
           lnB0=this.priorParameters(1);
           lnBvar=this.priorParameters(2)^2;           
           lnN0=this.priorParameters(3);
           lnNvar=this.priorParameters(4)^2;           
           lnS0=this.priorParameters(5);
           lnSvar=this.priorParameters(6)^2;           
           y = - 1/2 *(lnB - lnB0).^2./lnBvar-0.5*log(2*pi*lnBvar)...
               - 1/2 *(lnN - lnN0).^2./lnNvar-0.5*log(2*pi*lnNvar)...
               - 1/2 *((lnS1 - lnS0).^2+(lnS2 - lnS0).^2)./lnSvar...
               -log(2*pi*lnSvar);

           dy=zeros(size(param));
           dy(3)=- (lnB - lnB0)./lnBvar;
           dy(4)=- (lnN - lnN0)./lnNvar;
           dy(5)=- (lnS1 - lnS0)./lnSvar;
           dy(6)=- (lnS2 - lnS0)./lnSvar;
        end
        
        % check that the prior parameters are correct size
        function flag = hasValidPrior(this)
            flag=true;
            if(numel(this.priorParameters) ~=2 )
                flag=false;
            end          
        end
    end    
end