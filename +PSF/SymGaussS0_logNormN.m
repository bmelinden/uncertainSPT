classdef SymGaussS0_logNormN < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0), with
    % log-normal amplitude prior, i.e., normal prior on the fit parameter
    % lnN, with mean value parameters lnN0 and std-parameters lnNstd. 
    %
    % Construction:
    % P = SymGaussS0_logNormN('initialGuess',[mux muy lnB lnN lnS],'S0','priorParameters',[lnN0 lnNstd]),
    % P = SymGaussS0_logNormN('initialGuess',[mux muy lnB lnN lnS],'lambda',lambda,'NA',NA,'priorParameters',[lnN0 lnNstd]),
    % (the initial guess, NA,lambda or S0 are passed on to the SymGaussS0 constructor)
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = SymGaussS0_logNormN(varargin)
            % call superclass constructor
            this@PSF.SymGaussS0(varargin{:});
            % sanity check on priors
            if(numel(this.priorParameters)~=2 )
               error('PSF.SymGaussS0_logNormB needs 2 prior parameters [lnN0 lnNstd]')
            end
        end
        % log prior density
        function [y,dy] = logPrior(~,param)
            
            lnN=param(4);
            lnN0=this.priorParameters(1);
            lnNvar=this.priorParameters(2)^2;
           
            y= -1/2 *(lnN - lnN0).^2./lnNvar -0.5*log(2*pi*lnNvar);
           
            dy=zeros(size(param));
            dy(4)=- (lnN - lnN0)./lnNvar;
        end
    end
end