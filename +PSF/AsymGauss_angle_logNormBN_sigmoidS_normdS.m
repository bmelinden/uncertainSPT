classdef AsymGauss_angle_logNormBN_sigmoidS_normdS < PSF.AsymGauss_angle
    % An asymmetric Gaussian PSF model (see AsymGauss_angle) where
    %
    % * background B and amplitude N have log-normal priors with location
    % parameters lnB0, lnN0 and scale parameters lnBstd, lnNstd. Scale
    % parameter=inf means a flat prior on that variable.
    %
    % * the absolute difference dS=S1-S2=exp(lnS1)-exp(lnS2) has zero-mean
    % normal prior with standard deviation dSstd.
    % * the geometric mean PSF width, Sm=sqrt(S1*S2)=exp(lnSm), has
    % an unnormalized sigmoidal (soft step) prior at S0, with sharpness
    % dSscale, so that 
    % p(Sm)  = 1/(1+exp(-(Sm/S0-1)/dsScale))/Sm, or
    % p(lnSm)= 1/(1+exp(-(exp(lnSm)/S0-1)/dsScale))
    % This effectively enforces a lower bound on Sm, which can be chosen to
    % be that of a well-focused spot, S0 = 0.21*lambda/NA, and so the prior
    % can be constructed in two possible ways.
    % The factor 1/S is needed in order that the prior on log(Sm) must not
    % grow as Sm>>S0, which could potentially cause lnSm to grow without
    % bound during the optimization. dsScale=inf gives a constant prior
    %
    % Since the fit parameters are lnS1,lnS2, the prior picks up a factor
    % Sm*dS, so that
    % p0(lnS1,lnS2) = Sm*dS/2*N(dS;dSStd)*p0(Sm) 
    %               = dS/2*N(dS(lnS1,lnS2);dSStd)*p0(lnSm(lnS1,lnS2))
    %
    % Construction:
    % P = AsymGauss_angle_logNormBNdS_sigmoidS(...
    %  'priorParameters',pp, 'initialGuess',[mux muy lnB lnN lnS1 lnS2 v]),
    %  where the prior parameters can be either
    % pp = [lnB0 lnBstd lnN0 lnNstd dSstd dsScale S0], or
    % pp = [lnB0 lnBstd lnN0 lnNstd dSstd dsScale lambda NA], in
    % which case S0 = 0.21*lambda/NA is computed internally.
    
    properties (Constant)
        priorName='log-normal parameters';
    end
    properties
        priorParameters	=[];
        lambda=[];
        NA=[];
    end
    
    methods
        % Constructor
        function this = AsymGauss_angle_logNormBN_sigmoidS_normdS(varargin)
            % call superclass constructor
            this@PSF.AsymGauss_angle(varargin{:});
            % sanity check on prior parameters
            if(numel(this.priorParameters)==7) % the S0 is already given
            elseif(numel(this.priorParameters)==8) % then lambda,NA is given
                this.lambda=this.priorParameters(7);
                this.NA    =this.priorParameters(8);
                
                S0=0.21*this.lambda/this.NA;
                this.priorParameters=[this.priorParameters(1:6) S0];
            else
                error('PSF.AsymGauss_angle_logNormBNSR_sigmoidSM needs 7 or 8 prior parameters')
            end
            if( ~isempty(find(this.priorParameters([2 4 6]) <=0,1)) )
                error('PSF.AsymGauss_angle_logNormBNSR_sigmoidSM needs positive scale and std parameters.')
            end
        end
        
        % log prior distribution
        function [y,dy] = logPrior(this,param)
           y = 0;
           dy=zeros(size(param));
           % pp = [lnB0 lnBstd lnN0 lnNstd dSstd dsScale S0], or
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
 
           % [mux muy lnB lnN lnS1 lnS2 v]
           lnS1=param(5);
           lnS2=param(6);
           dSvar=this.priorParameters(5)^2;           
           if(isfinite(dSvar))               
               S1=exp(lnS1);
               S2=exp(lnS2);
               log_abs_dS=max(lnS1,lnS2)+log1p(-exp(-abs(lnS1-lnS2))); % log(|dS|) without cancellation               
               y      =y + log_abs_dS-log(2)-0.5*log(2*pi*dSvar) -1/2 *(S1-S2).^2/dSvar;
               
               ddS_dlnS1=S1./(S1-S2); % d|dS|/ dlnS1
               ddS_dlnS2=S2./(S2-S1); % d|dS|/ dlnS2
               
               dy(5)  = ddS_dlnS1 - S1.*(S1-S2)/dSvar;
               dy(6)  = ddS_dlnS2 - S2.*(S2-S1)/dSvar;
           end
           
           dsScale=this.priorParameters(6);
           if(isfinite(dsScale))
               S0=this.priorParameters(7);
               lnSm=(lnS1+lnS2)/2;
               Sm=exp(lnSm);
            
               y=-log1p(exp(-(Sm/S0-1)/dsScale));
               
               dydlnSm=Sm/dsScale/S0./(1+exp((Sm/S0-1)/dsScale));
               dy(5)  =dy(5)  +dydlnSm/2;
               dy(6)  =dy(6)  +dydlnSm/2;
           end
        end
    end    
end
