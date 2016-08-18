classdef SymGauss_logNormB < SymGauss
    
    properties
        prior = [];
    end
    
    methods
        % Constructor
        function this = SymGauss_logNormB(varargin)
            this.priorName = 'logNormB';
            
        end
        
        function [y] = logPrior(param)
           [~, ~, lnB] = translateFitParameters(param);
           
           lnB0=this.priorParameters(1);
           lnBstd=this.priorParameters(2);
           
           y = - 1/2 *(lnB - lnB0).^2./(lnBstd^2);
        end
        
    end
    
end