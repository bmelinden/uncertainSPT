classdef SymGauss_MLE < PSF.SymGauss
    % a symmetric Gaussian PSF model (SymGauss) with no prior
    %
    % Construction:
    % P = SymGauss_MLE('initialGuess',[mux muy lnB lnN lnS]),
    % P = SymGauss_MLE('initialGuess',[mux muy lnB lnN lnS],'priorParameters',[]),
    % (the initial guess is passed on to the SymGauss constructor)
    properties (Constant)
        priorName='no prior';
    end
    properties
        priorParameters	=[];
    end
    
    methods
        % Constructor
        function this = SymGauss_MLE(varargin)
            % call superclass constructor
            this@PSF.SymGauss(varargin{:});
            % sanity check on priors
            if(numel(this.priorParameters)>0)
               error('PSF.SymGauss_MLE needs empty or no priorParameters')
            end
        end
        % log prior density
        function [y,dy] = logPrior(~,param)
           y=0;
           dy=zeros(size(param));
        end
    end
end