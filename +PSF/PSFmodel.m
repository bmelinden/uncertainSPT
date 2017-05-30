classdef (Abstract) PSFmodel 
    % PSFmodel_subclass('initialGuess',pInit,'priorParameters',p0param)
    %
    % Abstract superclass for PSF model objects. All subclasses are to be 
    % initialized with the above syntax (and simply pass the constructor
    % arguments on upwards), but the density function, prior distribution,
    % and some helper functions need to be defined in supclasses.
    % Suggested structure: subclasses define the PSF model and parameters,
    % and different subclasses of these define prior distributions.
    properties (Abstract, Constant)
        modelName;
        priorName;
    end
    properties (Abstract)
        initialGuess;
        priorParameters;
    end
    methods 
        function this=PSFmodel(varargin)
            for n=1:2:nargin
                this.(varargin{n})=varargin{n+1};
            end
        end
    end
    methods (Abstract, Access = public)
        % convert fit parameters to a struct form: this should be
        % overloaded by the subclass defining the PSF model
        [outStruct,outPar] = convertToOutStruct(this,inPar, inStruct)
        [E,gradE] = psfDensity(this,x,y,param)
        [lnL0,dlnL0dp]= logPrior(this,param)  
    end
end

% overwrite +EMCCDfit/logL_psf.m with ClassDef/logL_psf2.m and update
% syntax in various places
% ML: make sure one can still recreate a pixellated image from fit
% parameters

% example using refineSingleFrame

% ML thought 1: should perhaps the convertToOutStruct method also do some
% unit conversion, for example rescale the length unit to ,e.g., nm? 
% In that case, a scaling factor has to be supplied somehow. I would
% suggest as a construction argument, since the user presumably knows what
% units fitting will be done in, and what units one wants the final answer
% to be.

% ML though 2: I did not set up access functions (e.g.
% this.getInitialGuess(), etc), but have no strong feelings about it.

% to use gradients: one will need to compute derivatives of the likelihood
% too in order to use this. This is not so easy to to efficiently, but it
% is still good to have the partial derivatives in place, given that it is
% a bit of a pain to dervive and validate them. And coputing gradients
% would speed things up quite a bit...
