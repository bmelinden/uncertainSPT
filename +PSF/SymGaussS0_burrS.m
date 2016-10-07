classdef SymGaussS0_burrS < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0), and 
    % Burr distribution prior (see documentation of the matlab class
    % prob.BurrDistribution) on the rescaled PSF width
    % y=(S-S0)/S0=exp(lnS), which means that the prior on lnS is  
    %
    % p(lnS) = pB.pdf(exp(lnS))*exp(lnS),
    %
    % where pB is a prob.BurrDistribution object with parameters alpha,c,k.
    %
    % Construction:
    % P = SymGaussS0_burrS('initialGuess',[mux muy lnB lnN lnS],'S0','priorParameters',[alpha c k]),
    % P = SymGaussS0_burrS('initialGuess',[mux muy lnB lnN lnS],'lambda',lambda,'NA',NA,'priorParameters',[alpha c k]),
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
        function this = SymGaussS0_burrS(varargin)
            % call superclass constructor
            this@PSF.SymGaussS0(varargin{:});
            % sanity check on priors
            if(numel(this.priorParameters)~=3 )
               error('PSF.SymGaussS0_burrS needs 3 prior parameters [alpha c k]')
            end
            this.pB=prob.BurrDistribution(this.priorParameters(3),this.priorParameters(4),this.priorParameters(5));
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
            lnS=param(5);
            [lnPB,dlnPB]=lnP_expBurr(lnS);
            
            y=lnPB;
           
            dy=zeros(size(param));
            dy(5)= dlnPB;
        end
    end
end
