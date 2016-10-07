classdef SymGaussS0_logNormBN < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0), with
    % log-normal background amplitude priors, i.e., normal priors on the
    % fit parameters lnB, lnN, with mean value parameters lnN0, lnB0 and
    % std-parameters lnNstd, lnBstd.
    %
    % Construction:
    % P = SymGaussS0_logNormBN('initialGuess',[mux muy lnB lnN lnS],'S0','priorParameters',[lnB0 lnBstd lnN0 lnNstd]),
    % P = SymGaussS0_logNormBN('initialGuess',[mux muy lnB lnN lnS],'lambda',lambda,'NA',NA,'priorParameters',[lnB0 lnBstd lnN0 lnNstd]),
    % (the initial guess, NA,lambda or S0 are passed on to the SymGaussS0 constructor)
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = SymGaussS0_logNormBN(varargin)
            % call superclass constructor
            this@PSF.SymGaussS0(varargin{:});
            % sanity check on priors
            if(numel(this.priorParameters)~=4 )
               error('PSF.SymGaussS0_logNormBN needs 4 prior parameters [lnB0 lnBstd lnN0 lnNstd]')
            end
        end
        % log prior density
        function [y,dy] = logPrior(~,param)
            
            lnB=param(3);
            lnB0=this.priorParameters(1);
            lnBvar=this.priorParameters(2)^2;

            lnN=param(4);
            lnN0=this.priorParameters(3);
            lnNvar=this.priorParameters(4)^2;
           
            y= -1/2 *(lnB - lnB0).^2./lnBvar -0.5*log(2*pi*lnBvar)...
               -1/2 *(lnN - lnN0).^2./lnNvar -0.5*log(2*pi*lnNvar);
           
            dy=zeros(size(param));
            dy(3)=- (lnB - lnB0)./lnBvar;
            dy(4)=- (lnN - lnN0)./lnNvar;
        end
    end
end