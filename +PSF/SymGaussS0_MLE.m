classdef SymGaussS0_MLE < PSF.SymGaussS0
    % a symmetric Gaussian PSF model with width offset (SymGaussS0), no prior
    %
    % Construction:
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lnS],'Smin',Smin),
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lnS],'lambda',lambda,'NA',NA),
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lnS],'Smin','priorParameters',[]),
    % P = SymGaussS0_MLE('initialGuess',[mux muy lnB lnN lnS],'lambda',lambda,'NA',NA,'priorParameters',[]),
    % (the initial guess, NA,lambda,Smin are passed on to the SymGaussS0 constructor)
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