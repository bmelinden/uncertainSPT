classdef AsymGauss_lndS_angle_BNlnN_dSexp_N < PSF.AsymGauss_lndS_angle
    % An asymmetric Gaussian PSF model with PSF width offset (see
    % PSF.AsymGauss_lndS_angle) parameterized by fit parameters
    % lnB,lnN,lndS1,lndS2,v.
    % 
    % lnB, lnN: logNormal with scale parameters lnB0, lnN0, and shape
    % parameters lnBstd, lnNstd. shape parameter = inf means no prior 
    %
    % lndS1,lndS2: 
    % 1) (dS1+dS2)/2 ~ exp(S0), exp-distributed with mean value dS0, and 
    % 2) S1-S2   ~ N(0,dSstd), i.e.
    %    p0(lndS1,lndS2) =
    %    1/dS0*exp(lndS1 + lndS2 -(exp(lndS1)+exp(lndS2))/2/dS0)
    %   *1/sqrt(2*pi*dSstd^2)*exp(-0.5*(exp(lndS1)-exp(lndS2))^2/dSstd^2)
    % S0, dSstd = inf means the corresponding contribution is flat (with
    % the exp(lndS1+lndS2)-factor "belonging" to the exp-distribution) and
    % if S0=dSstd=inf, then p0=const. Note that this prior is not properly
    % normalized, but this does not influence the shape of the posterior
    % distribution. It does mean that one cannot optimize prior parameters
    % however.
    %
    % Note: in the limit dSstd -> 0, this prior approaches the symmetric
    % prior with an exponential prior dS with mean value in actual
    % length units (not in units of S0).
    %
    % Construction:
    % P=PSF.AsymGauss_lndS_angle_BNlnN_dSexp_N('lambda',lambda,'NA',NA,...
    %           'initialGuess',[mux muy lnB lnN lndS1 lndS2 v],...
    %           'priorParameters',[lnB0 lnBstd lnN0 lnNstd dS0 dSstd])
    % P=PSF.AsymGauss_angle_MLE('S0',S0,...
    % (The initial guess is passed on to the AsymGauss_lndS_angle constructor.)
    properties (Constant)
        priorName='N,B~lnN, dS1,dS2~exp*N';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = AsymGauss_lndS_angle_BNlnN_dSexp_N(varargin)
            % call superclass constructor
            this@PSF.AsymGauss_lndS_angle(varargin{:});

            % sanity check on priors
            if(numel(this.priorParameters)~=6)
               error('PSF.AsymGauss_lndS_angle_BNlnN_dSexp_N needs 6 prior parameters')
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
            dS1=exp(lndS1);
            dS2=exp(lndS2);
            dS0=this.priorParameters(5);
            dSvar=this.priorParameters(6);
            if(~isfinite(dS0) && ~isfinite(dSvar))
                % nothing to do
            else
                %    p0(lndS1,lndS2) =
                %    1/dS0*exp(lndS1 + lndS2 -(exp(lndS1)+exp(lndS2))/2/dS0)
                %   *1/sqrt(2*pi*dSstd^2)*exp(-0.5*(exp(lndS1)-exp(lndS2))^2/dSstd^2)
                if(isfinite(dS0))
                    y=y+lndS1+lndS2-log(dS0)-(dS1+dS2)/2/dS0;
                    dy(5)=1+dy(5)-dS1/dS0/2;
                    dy(6)=1+dy(6)-dS2/dS0/2;
                end
                if(isfinite(dSvar))
                    y=y-0.5*log(2*pi*dSvar)-0.5*(dS1-dS2)^2/dSvar;
                    dy(5)=dy(5)-(dS1-dS2)*dS1/dSvar;
                    dy(6)=dy(6)-(dS2-dS1)*dS2/dSvar;
                end
            end
        end
    end
end
