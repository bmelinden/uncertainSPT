classdef SymGaussS0_logNormBN_burrS < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0), with
    % log-normal background amplitude priors, i.e., normal priors on the
    % fit parameters lnB, lnN, with mean value parameters lnN0, lnB0 and
    % std-parameters lnNstd, lnBstd, and Burr distribution (see
    % documentation of the matlab class prob.BurrDistribution) on the
    % rescaled PSF width y=(S-S0)/S0=exp(lnS), which means that the prior
    % on lnS is
    %
    % p(lnS) = pB.pdf(exp(lnS))*exp(lnS),
    %
    % where pB is a prob.BurrDistribution object with parameters alpha,c,k.
    %
    % Construction:
    % P = SymGaussS0_logNormBN_burrS('initialGuess',[mux muy lnB lnN lnS],'S0','priorParameters',[lnB0 lnBstd lnN0 lnNstd alpha c k]),
    % P = SymGaussS0_logNormBN_burrS('initialGuess',[mux muy lnB lnN lnS],'lambda',lambda,'NA',NA,'priorParameters',[lnB0 lnBstd lnN0 lnNstd  alpha c k]),
    % (the initial guess, NA,lambda or S0 are passed on to the SymGaussS0 constructor)
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
        pB=[];
    end
    
    methods
        % Constructor
        function this = SymGaussS0_logNormBN_burrS(varargin)
            % call superclass constructor
            this@PSF.SymGaussS0(varargin{:});
            % sanity check on priors
            if(numel(this.priorParameters)~=7 )
               error('PSF.SymGaussS0_logNormBN_burrS needs 7 prior parameters [lnB0 lnBstd lnN0 lnNstd alpha c k]')
            end
            this.pB=prob.BurrDistribution(this.priorParameters(5),this.priorParameters(6),this.priorParameters(7));
        end
        % Burr distribution log likelihood
        function [lnP,dlnPdlnS] = lnP_expBurr(this,lnS)
            
            a=this.pB.alpha;
            k=this.pB.k;
            c=this.pB.c;
            
            lnP = lnS + log(k*c/a) ...
                +(c-1)*(lnS-log(a)) ...
                -(k+1)*log1p((exp(lnS)/a).^c);
            
            dlnPdlnS = c -(k+1)./(1+exp(c*lnS)/a^c)*c/a^c.*exp(c*lnS);
        end
        % log prior density
        function [y,dy] = logPrior(~,param)
            
            lnB=param(3);
            lnB0=this.priorParameters(1);
            lnBvar=this.priorParameters(2)^2;

            lnN=param(4);
            lnN0=this.priorParameters(3);
            lnNvar=this.priorParameters(4)^2;
            
            lnS=param(5);
            [lnPB,dlnPB]=lnP_expBurr(lnS);
            
            y= -1/2 *(lnB - lnB0).^2./lnBvar -0.5*log(2*pi*lnBvar)...
               -1/2 *(lnN - lnN0).^2./lnNvar -0.5*log(2*pi*lnNvar)...
               +lnPB;
           
            dy=zeros(size(param));
            dy(3)=- (lnB - lnB0)./lnBvar;
            dy(4)=- (lnN - lnN0)./lnNvar;
            dy(5)= dlnPB;
        end
    end
end