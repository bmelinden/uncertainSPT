classdef AsymGauss_angle_MLE < AsymGauss_angle
    % An asymmetric Gaussian PSF model (see AsymGauss_angle) where the
    % background  intensity B has a log-normal prior (or log(B) is Normal)
    % with location parameter lnB0 and scale parameter lnBstd, which means
    % that <ln B> = lnB0, Std(ln B) lnBstd, <B> = exp(lnB0+lnBstd^2/2),
    % Std(B) = exp(lnB0+lnBstd^2/2)*sqrt(exp(lnBstd^2)-1),
    % see https://en.wikipedia.org/wiki/Log-normal_distribution.
    %
    % Construction:
    % P = AsymGauss_angle_MLE('initialGuess',[mux muy lnB lnN lnS1 lnS2 v]
    % (The initial guess is passed on to the AsymGauss_angle constructor.)
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = AsymGauss_angle_MLE(varargin)
            % call superclass constructor
            this@AsymGauss_angle(varargin{:});

            % check validity of the model
            if( ~this.hasValidInitialGuess() || ~this.hasValidPrior())
                disp('AsymGauss_angle_MLE constructor args:')
                varargin{:} %#ok<NOPRT>
                this %#ok<NOPRT>
                error('AsymGauss_angle_MLE not correctly set up.')
            end
        end
        
        % log prior distribution
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