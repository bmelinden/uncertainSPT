classdef AsymGaussS0_BNlnN_dSexp_N < PSF.AsymGaussS0
    % An asymmetric Gaussian PSF model with PSF width offset (see
    % PSF.AsymGaussS0) parameterized by fit parameters
    % lnB,lnN,lndS1,lndS2,v.
    % 
    % lnB, lnN: logNormal with scale parameters lnB0, lnN0, and shape
    % parameters lnBstd, lnNstd. shape parameter = inf means no prior 
    %
    % lndS1,lndS2: 
    % dS1=exp(lndS1)=(S1-S0)/S0
    % dS2=exp(lndS2)=(S2-S0)/S0    
    % 1) (dS1+dS2)/2 ~ exp(1/dS0), exp-distributed with mean value dS0, and 
    % 2) ln(S2/S1)   ~ N(0,dSstd), i.e.
    %    p0(lndS1,lndS2) =
    %        [ exp(lndS1)/(1+exp(-lndS2)) + exp(lndS2)/(1+exp(-lndS1)) ]
    %       *1/dS0*exp(-(exp(lndS1)+exp(lndS2))/2/dS0)
    %       *1/sqrt(2*pi)/dSstd*exp(...
    %           -0.5*( ( ln(1+exp(lndS2))-ln(1+exp(lndS1)) )/dSstd )^2 )
    % S0, dSstd must be finite.
    %
    % Construction:
    % P=PSF.AsymGaussS0_BNlnN_dSexp_N('lambda',lambda,'NA',NA,...
    %           'initialGuess',[mux muy lnB lnN lndS1 lndS2 v],...
    %           'priorParameters',[lnB0 lnBstd lnN0 lnNstd dS0 dSstd])
    % or
    % P=PSF.AsymGaussS0_BNlnN_dSexp_N('S0',S0,...
    % (The initial guess is passed on to the AsymGaussS0 constructor.)
    properties (Constant)
        priorName='N,B~lnN, dS1,dS2~exp*N';
    end
    properties
        priorParameters	=[];
    end    
    methods
        % Constructor
        function this = AsymGaussS0_BNlnN_dSexp_N(varargin)
            % call superclass constructor
            this@PSF.AsymGaussS0(varargin{:});
            
            % sanity check on priors
            if(numel(this.priorParameters)~=6)
                error('PSF.AsymGaussS0_BNlnN_dSexp_N needs 6 prior parameters')
            end
            if( ~isfinite(this.priorParameters(5)) || ~isfinite(this.priorParameters(6)) )
                error('PSF.AsymGaussS0_BNlnN_dSexp_N : prior parameters dS0 dSstd must be finite.')
            end
        end
        
        % log prior distribution
        function [y,dy] = logPrior(this,param)
            
            % ML: verified numerically 2016-11-08
            y = 0;
            dy=zeros(size(param));
            % pp = [lnB0 lnBstd lnN0 lnNstd dS0 dSstd]
            lnBvar =this.priorParameters(2)^2;
            if(isfinite(lnBvar))
                lnB    =param(3);
                lnB0   =this.priorParameters(1);
                y      =y-1/2*(lnB - lnB0).^2./lnBvar -0.5*log(2*pi*lnBvar);
                dy(3)  =-(lnB - lnB0)./lnBvar;
            end
            
            lnNvar =this.priorParameters(4)^2;
            if(isfinite(lnNvar))
                lnN    =param(4);
                lnN0   =this.priorParameters(3);
                y      =y-1/2 *(lnN - lnN0).^2./lnNvar -0.5*log(2*pi*lnNvar);
                dy(4)  =-(lnN - lnN0)./lnNvar;
            end
            
            lndS1=param(5);
            lndS2=param(6);
            dS0=this.priorParameters(5);
            dSvar=this.priorParameters(6)^2;
            
            dS1=exp(lndS1);
            dS2=exp(lndS2);
            
            lns1=log1p(dS1);
            lns2=log1p(dS2);
            
            %    p0(lndS1,lndS2) =
            %        [ exp(lndS1)/(1+exp(-lndS2)) + exp(lndS2)/(1+exp(-lndS1)) ]
            %       *1/dS0*exp(-(exp(lndS1)+exp(lndS2))/2/dS0)
            %       *1/sqrt(2*pi)/dSstd*exp(...
            %           -0.5*( ( ln(1+exp(lndS2))-ln(1+exp(lndS1)) )/dSstd )^2 )
            y=y-log(2)+lndS1+lndS2+log(1/(1+dS1)+1/(1+dS2))...
                -log(dS0)-(dS1+dS2)/2/dS0 ...
                -0.5*log(2*pi*dSvar)-0.5*(lns1-lns2)^2/dSvar;
            
            % partial derivatives verified numerically 2016-11-21/ML
            dy(5)=dy(5)+1-dS1/(1+dS1+(1+dS1)^2/(1+dS2)); % d ln|dz/dy| / dy1
            dy(6)=dy(6)+1-dS2/(1+dS2+(1+dS2)^2/(1+dS1));
            
            dy(5)=dy(5)-dS1/dS0/2; % exp-contributions
            dy(6)=dy(6)-dS2/dS0/2;
            
            dy(5)=dy(5)-(lns1-lns2)/dSvar*dS1/(1+dS1); % logNorm-contributions
            dy(6)=dy(6)-(lns2-lns1)/dSvar*dS2/(1+dS2);
        end
    end
end
