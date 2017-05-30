classdef SymGaussS0_MLE < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0) and
    % excess psf width dS=(S-S0)/S0=exp(lndS).
    %
    % Construction:
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lndS],'Smin',Smin),
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lndS],'lambda',lambda,'NA',NA),
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lndS],'Smin','priorParameters',[]),
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lndS],'lambda',lambda,'NA',NA,'priorParameters',[]),
    % (the initial guess, NA,lambda,Smin are passed on to the SymGaussS0 constructor)
    %
    % Since fitting is done using the lndS parameter, there can be numerical
    % difficulties to get good Lapace error estimates using this PSF model
    % when lndS<<0. In this case, it is probably better to use a soft
    % cut-off instead, e.g., SymGauss_sigmoidS.m
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = SymGaussS0_MLE(varargin)
            % call superclass constructor
            this@PSF.SymGaussS0(varargin{:});
            % sanity check on priors
            if(numel(this.priorParameters)>0)
               error('PSF.SymGaussS0_MLE needs empty or no priorParameters')
            end
        end
        % log prior density
        function [y,dy] = logPrior(~,param)
           y=0;
           dy=zeros(size(param));
        end
    end
end
