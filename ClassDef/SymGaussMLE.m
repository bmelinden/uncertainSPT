classdef SymGaussMLE < SymGauss
    
    properties
        prior = [];
    end
    
    methods
        % Constructor
        function this = SymGaussMLE(varargin)
            this.priorName = 'MLE';
        end
        
        function y = logPrior(this, param)
           y=0;
        end
    end
    
end