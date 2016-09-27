classdef flatBG_MLE < PSF.flatBG
    % A flat background model (flatBG), where the log background intensity
    % lnB has a flat prior.
    %
    % Construction:
    % P = flatBG_logNormB('initialGuess',lnB,'priorParameters',[]),
    % (the initial guess is passed on to the flatBG constructor)
    properties (Constant)
        priorName='lnN-background';
    end
    properties
        priorParameters	=[];
    end
    methods
        % Constructor
        function this = flatBG_MLE(varargin)
            % call superclass constructor
            this@PSF.flatBG(varargin{:});
            % sanity check on priorParameters
            if(~isempty(this.priorParameters)>0)
                error('PSF.flatBG_MLE accepts no priorParameters.')
            end
        end
        % compute log prior density
        function [y,dy] = logPrior(this,param)
           y = 0;
           dy= 0;
        end        
    end
end