classdef AsymGauss_angle_logNormBNSR_sigmoidSM < PSF.AsymGauss_angle
    % An asymmetric Gaussian PSF model (see AsymGauss_angle) where
    %
    % * background B, amplitude N, and principal width ration SR=S1/S2 have
    % log-normal priors with location parameters lnB0, lnN0, lnSR0=0, and
    % scale parameters lnBstd, lnNstd, lnSRstd. Scale parameter=inf means
    % a flat prior on that variable.
    %
    % * the geometric mean PSF width, Sm=sqrt(S1*S2)=exp(lnSm), has
    % an unnormalized sigmoidal (soft step) prior at S0, with sharpness
    % dSscale, so that 
    % p(Sm)  = 1/(1+exp(-(Sm/S0-1)/dsScale))/S, or
    % p(lnSm)= 1/(1+exp(-(exp(lnSm)/S0-1)/dsScale))
    % This effectively enforces a lower bound on Sm, which can be chosen to
    % be that of a well-focused spot, S0 = 0.21*lambda/NA, and so the prior
    % can be constructed in two possible ways.
    % The factor 1/S is needed in order that the prior on log(Sm) must not
    % grow as Sm>>S0, which could potentially cause lnSm to grow without
    % bound during the optimization. dsScale=inf gives a constant prior
    %
    % Construction:
    % P = AsymGauss_angle_logNormBNdS_sigmoidS(...
    %  'priorParameters',pp, 'initialGuess',[mux muy lnB lnN lnS1 lnS2 v]),
    %  where the prior parameters can be either
    % pp = [lnB0 lnBstd lnN0 lnNstd lnSRstd dsScale S0], or
    % pp = [lnB0 lnBstd lnN0 lnNstd lnSRstd dsScale lambda NA], in
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
        function this = AsymGauss_angle_logNormBNSR_sigmoidSM(varargin)
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
                error('PSF.AsymGauss_angle_logNormBNdS_sigmoidS needs positive scale parameters.')
            end
        end
        
        % log prior distribution
        function [y,dy] = logPrior(this,param)
           y = 0;
           dy=zeros(size(param));
           % pp = [lnB0 lnBstd lnN0 lnNstd lnSRstd dsScale S0], or
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
           lnSRvar=this.priorParameters(5)^2;
           if(isfinite(lnSRvar))               
               y      =y-1/2 *(lnS1 - lnS2).^2/lnSRvar -0.5*log(2*pi*lnSRvar);
               dy(5)  = -(lnS1 - lnS2)/lnSRvar;
               dy(6)  = -(lnS2 - lnS1)/lnSRvar;
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
