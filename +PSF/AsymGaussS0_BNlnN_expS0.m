classdef AsymGaussS0_BNlnN_expS0 < PSF.AsymGaussS0
    % An asymmetric Gaussian PSF model with PSF width offset (see
    % PSF.AsymGaussS0) parameterized by fit parameters
    % lnB,lnN,lndS1,lndS2,v.
    % 
    % lnB, lnN: logNormal with scale parameters lnB0, lnN0, and shape
    % parameters lnBstd, lnNstd. shape parameter = inf means no prior 
    %
    % lndS1,lndS2: exponential priors on the excess psf widths
    % dSj=(Sj-S0)/S0=exp(lndSj), j=1,2, which translates to a prior density  
    %
    % p0(lndSj) = 1/dS0*exp(lndSj-exp(lndSj)/dS0)
    % 
    % for an exponential prior with mean value dS0. dS0=inf give a flat
    % prior on lndSj, i.e., p0(dSj)~1/dSj.
    %
    % Construction:
    % P=PSF.AsymGaussS0_BNlnN_dSexp_N('lambda',lambda,'NA',NA,...
    %           'initialGuess',[mux muy lnB lnN lndS1 lndS2 v],...
    %           'priorParameters',[lnB0 lnBstd lnN0 lnNstd dS0])
    % or
    % P=PSF.AsymGaussS0_BNlnN_dSexp_N('S0',S0,...
    % (The initial guess is passed on to the AsymGaussS0 constructor.)
    properties (Constant)
        priorName='N,B~lnN, dSj~exp';
    end
    properties
        priorParameters	=[];
    end    
    methods
        % Constructor
        function this = AsymGaussS0_BNlnN_expS0(varargin)
            % call superclass constructor
            this@PSF.AsymGaussS0(varargin{:});
            
            % sanity check on priors
            if(numel(this.priorParameters)~=5)
                error('PSF.AsymGaussS0_BNlnN_expS0 needs 5 prior parameters')
            end
        end
        
        % log prior distribution
        function [y,dy] = logPrior(this,param)
            
            % ML: verified numerically 2016-11-08
            y = 0;
            dy=zeros(size(param));
            % pp = [lnB0 lnBstd lnN0 lnNstd dS0 SrStd]
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
            if(isfinite(dS0))
                y=y-log(dS0)+lndS1-exp(lndS1)/dS0;
                dy(5) = 1-exp(lndS1)/dS0;
                y=y-log(dS0)+lndS2-exp(lndS2)/dS0;
                dy(6) = 1-exp(lndS2)/dS0;
            end
        end
    end
end
