classdef SymGaussS0_logNormBN_burrS < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0),
    % parameterized by the log-excess PSF width lndS, exp(lndS)=dS=(S-S0)/S0,
    % where S0 is the minimum PSF width.
    %
    % Priors:
    % - log-normal background amplitude priors, i.e., normal priors on the
    % fit parameters lnB, lnN, with mean value parameters lnN0, lnB0 and
    % std-parameters lnNstd, lnBstd. Setting
    % lnBstd, lnNstd to inf gives flat priors on lnB or lnN, respectively.
    %
	% - Burr distribution (see documentation of the matlab class
	% prob.BurrDistribution) on the excess PSF width dS=(S-S0)/S0=exp(lndS),
	% which means that the prior on lndS is
    %
    % p(lndS) = pB.pdf(exp(lndS))*exp(lndS),
    %
    % where pB is a prob.BurrDistribution object with parameters alpha,c,k.
    %
    % Construction:
    % P = SymGaussS0_logNormBN_burrS('initialGuess',[mux muy lnB lnN lndS],'S0','priorParameters',[lnB0 lnBstd lnN0 lnNstd alpha c k]),
    % P = SymGaussS0_logNormBN_burrS('initialGuess',[mux muy lnB lnN lndS],'lambda',lambda,'NA',NA,'priorParameters',[lnB0 lnBstd lnN0 lnNstd  alpha c k]),
    % (the initial guess, NA,lambda or S0 are passed on to the SymGaussS0 constructor)
    properties (Constant)
        priorName='logNorm(B,N), burr(dS)';
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
        function [lnP,dlnPdlndS] = lnP_expBurr(this,lndS)
            
            a=this.pB.alpha;
            k=this.pB.k;
            c=this.pB.c;
            
            lnP = lndS + log(k*c/a) ...
                +(c-1)*(lndS-log(a)) ...
                -(k+1)*log1p((exp(lndS)/a).^c);
            
            dlnPdlndS = c -(k+1)./(1+exp(c*lndS)/a^c)*c/a^c.*exp(c*lndS);
        end
        % log prior density
        function [y,dy] = logPrior(this,param)
                        
            lndS=param(5);
            [lnPB,dlnPB]=this.lnP_expBurr(lndS);
            
            y=lnPB;
            dy=zeros(size(param));
            dy(5)= dlnPB;
           
            lnBvar=this.priorParameters(2)^2;
            if(isfinite(lnBvar))
                lnB=param(3);
                lnB0=this.priorParameters(1);
                y=y -1/2 *(lnB - lnB0).^2./lnBvar -0.5*log(2*pi*lnBvar);
                dy(3)=- (lnB - lnB0)./lnBvar;
            end
            
            lnNvar=this.priorParameters(4)^2;
            if(isfinite(lnNvar))
                lnN=param(4);
                lnN0=this.priorParameters(3);
                y=y -1/2 *(lnN - lnN0).^2./lnNvar -0.5*log(2*pi*lnNvar);
                dy(4)=- (lnN - lnN0)./lnNvar;
            end

            
            
        end
    end
end
