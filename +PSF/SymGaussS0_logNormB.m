classdef SymGaussS0_logNormB < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0), with
    % log-normal background prior, i.e., normal priors on the fit parameter
    % lnB, with mean value parameters lnB0 and std-parameters lnBstd. 
    %
    % Construction:
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lnS],'S0','priorParameters',[lnB0 lnBstd]),
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lnS],'lambda',lambda,'NA',NA,'priorParameters',[lnB0 lnBstd]),
    % (the initial guess, NA,lambda or S0 are passed on to the SymGaussS0 constructor)
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = SymGaussS0_logNormB(varargin)
            % call superclass constructor
            this@PSF.SymGaussS0(varargin{:});
            % sanity check on priors
            if(numel(this.priorParameters)~=2 )
               error('PSF.SymGaussS0_logNormB needs 2 prior parameters [lnB0 lnBstd]')
            end
        end
        % log prior density
        function [y,dy] = logPrior(~,param)
            
            lnB=param(3);
            lnB0=this.priorParameters(1);
            lnBvar=this.priorParameters(2)^2;
           
            y= -1/2 *(lnB - lnB0).^2./lnBvar -0.5*log(2*pi*lnBvar);
           
            dy=zeros(size(param));
            dy(3)=- (lnB - lnB0)./lnBvar;
        end
    end
end