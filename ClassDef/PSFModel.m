classdef (Abstract) PSFModel
    
    properties
        initialGuess    = [];
        priorParameters = [];
        priorName       = [];
    end
    
    methods
        % Constructor
        function this = PSFModel(varargin)
            for i = 1 : 2 : length(varargin)
                this.(varargin{i}) = varargin{i+1};
            end
        end
        
        function [outStruct] = convertToOutStruct(this, inStruct, inPar)
            
        end
    end
end