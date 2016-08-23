classdef SymGauss_MLE < SymGauss
    % a symmetric Gaussian PSF model (SymGauss) with no prior
    %
    % Construction:
    % P = SymGauss_MLE('initialGuess',[mux muy lnB lnN lnS]),
    % (the initial guess is passed on to the SymGauss constructor)
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = SymGauss_MLE(varargin)
            % call superclass constructor
            this@SymGauss(varargin{:});
            % check validity of the model
            if( ~this.hasValidInitialGuess() || ~this.hasValidPrior())
                disp('SymGauss_MLE constructor args:')
                varargin{:} %#ok<NOPRT>
                this %#ok<NOPRT>
                error('SymGauss_MLE not correctly set up.')
            end
        end
        
        % compute log prior density
        function [y,dy] = logPrior(this,param)
           y=0;
           dy=zeros(size(param));
        end
        
        % check that the prior parameters are correct size
        function flag = hasValidPrior(this)
            flag=true;                   
        end
        
    end
    
end