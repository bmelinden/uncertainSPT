classdef flatBG < PSF.PSFmodel
    % A flat background model, described by parameters
    % p = [lnB], and with density
    %
    % E = B, where B= exp(lnB),
    %
    % NOTE 1 : this model has no position parameters, and so cannot be used
    % with EMCCDfit.refineSingleFrame.
    % NOTE 2: this is a continuous model, so that exp(lnB) is the
    % background intensity per unit area, not per pixel. (But when doing
    % math in pixel units the area per pixel is 1, and then this does not
    % matter.)
    
    properties (Constant)
        modelName = 'Symmetric Gaussian';
    end
    properties
        initialGuess=0;
    end
    
    methods (Access = public)
        % Constructor
        function this = flatBG(varargin)
            this@PSF.PSFmodel(varargin{:});
            % sanity check for initialGuess
            if(numel(this.initialGuess)~=1)
                error('flatBG needs initialGuess scalar')
            end
        end
        % PSFmodel functions in this class level
        function [E,dEdp]= psfDensity(this, xx, yy, param)
            %   [E,dEdp]= psfDensity(this, xx, yy, param)
            E=exp(param(1))*ones(size(xx));
            dEdp=E;
        end
        
        function outStruct = convertToOutStruct(this,inPar, inStruct)
            % outStruct = convertToOutStruct(this,inPar, inStruct)
            % Convert psf fit parameters to a nice-looking parameter
            % struct
            %
            % inPar     : parameter vector to convert.
            % inStruct  : Optional struct to which parameter fields are
            %             added.
            %
            % output: parameter struct outStruct, with fields
            % outStruct.background  : background intensity photons/area = exp(lnB)
            
            if(nargin==3)
                outStruct = inStruct;
            else
                outStruct=struct;
            end
            outStruct.background=exp(inPar(1));
        end
    end
end
