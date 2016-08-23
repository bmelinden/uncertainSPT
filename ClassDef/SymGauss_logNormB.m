classdef SymGauss_logNormB < SymGauss
    % a symmetric Gaussian PSF model (SymGauss) where the background
    % intensity B has a log-normal prior with location parameter lnB0 and
    % scale parameter lnBstd, which means that
    % <ln B> = lnB0, Std(ln B) lnBstd, <B> = exp(lnB0+lnBstd^2/2),
    % Std(B) = exp(lnB0+lnBstd^2/2)*sqrt(exp(lnBstd^2)-1),
    % see https://en.wikipedia.org/wiki/Log-normal_distribution.
    %
    % Construction:
    % P = SymGauss_logNormB('priorParameters',[mux muy lnB lnN lnS],'initialGuess',[B0 N0 S0]),
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
            this@SymGauss(varargin{:});            
            % parse parametes for priorParameters
            k=0;
            while(k < nargin)
                k=k+1;
                vk=varargin{k};
                if(ischar(vk) && strcmp(vk,'priorParameters') )
                    % then store the prior parameters
                    k=k+1;
                    this.priorParameters = varargin{k};
                    constructor_satisfied=true;
                end
            end
            % check validity of the model
            if( ~this.hasValidInitialGuess() || ~this.hasValidPrior())
                disp('SymGauss_logNormB constructor args:')
                varargin{:} %#ok<NOPRT>
                this %#ok<NOPRT>
                error('SymGauss_logNormB not correctly set up.')
            end
        end
        
        % compute log prior density
        function [y,dy] = logPrior(this,param)
           [~, ~, lnB] = this.translateFitParameters(param);           
           lnB0=this.priorParameters(1);
           lnBstd=this.priorParameters(2);
           
           y = - 1/2 *(lnB - lnB0).^2./(lnBstd^2);
           dy=zeros(size(param));
           dy(3)=- (lnB - lnB0)./(lnBstd^2);
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