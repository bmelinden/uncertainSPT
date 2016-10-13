classdef SymGaussS0_burrS < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0),
    % parameterized by the log-excess PSF width lndS, exp(lndS)=dS=(S-S0)/S0,
    % where S0 is the minimum PSF width.
    %
    % Priors:
    % - Burr distribution (see documentation of the matlab class
    % prob.BurrDistribution) on the excess PSF width dS=(S-S0)/S0=exp(lndS),
    % which means that the prior on lndS is
    %
    % p(lndS) = pB.pdf(exp(lndS))*exp(lndS),
    %
    % where pB is a prob.BurrDistribution object with parameters alpha,c,k.
    %
    % Construction:
    % P = SymGaussS0_burrS('initialGuess',[mux muy lnB lnN lndS],'S0','priorParameters',[alpha c k]),
    % P = SymGaussS0_burrS('initialGuess',[mux muy lnB lnN lndS],'lambda',lambda,'NA',NA,'priorParameters',[alpha c k]),
    % (the initial guess, NA,lambda or S0 are passed on to the SymGaussS0 constructor)
    properties (Constant)
        priorName='burr(dS)';
    end
    properties
        priorParameters	=[];
        pB=[];
    end
    
    methods
        % Constructor
        function this = SymGaussS0_burrS(varargin)
            % call superclass constructor
            this@PSF.SymGaussS0(varargin{:});
            % sanity check on priors
            if(numel(this.priorParameters)~=3 )
                error('PSF.SymGaussS0_burrS needs 3 prior parameters [alpha c k]')
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
        end
    end
end
