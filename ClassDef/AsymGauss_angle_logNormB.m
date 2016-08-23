classdef AsymGauss_angle_logNormB < AsymGauss_angle
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
            this@AsymGauss_angle(varargin{:});
            % parse parametes for priorParameters
            k=0;
            while(k < nargin)
                k=k+1;
                vk=varargin{k};
                if(ischar(vk) && strcmp(vk,'priorParameters') )
                    % then store the prior parameters
                    k=k+1;
                    this.priorParameters = varargin{k};
                end
            end
            % check validity of the model
            if( ~this.hasValidInitialGuess() || ~this.hasValidPrior())
                disp('AsymGauss_angle_logNormB constructor args:')
                varargin{:} %#ok<NOPRT>
                this %#ok<NOPRT>
                error('AsymGauss_angle_logNormB not correctly set up.')
            end
        end
        
        % log prior distribution
        function [y] = logPrior(this,param)
           lnB=param(3);
           lnB0=this.priorParameters(1);
           lnBstd=this.priorParameters(2);           
           y = - 1/2 *(lnB - lnB0).^2./(lnBstd^2);
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