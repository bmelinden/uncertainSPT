classdef PSFModel
    % PSFmodel_subclass('initialGuess',pInit,'priorParameters',p0param)
    %
    % Abstract superclass for PSF model objects. All subclasses are
    % initialized with the above syntax (and simply pass the constructor
    % arguments on upwards), but the density function, prior distribution,
    % and some helper functions need to be defined in supclasses.
    % Suggested structure: subclasses define the PSF model and parameters,
    % and different subclasses of these define prior distributions.
    %
    properties (Abstract, Constant) % 
        modelName;
        priorName;
    end
    properties (Abstract)
        initialGuess;
        priorParameters;
    end
    methods (Abstract, Access = public)
        % convert fit parameters to a struct form: this should be
        % overloaded by the subclass defining the PSF model
        outStruct = convertToOutStruct(this,inPar, inStruct)
        [E,gradE] = psfDensity(this,x,y,param)
        [lnL0,dlnL0dp]= logPrior(this,param)  
        % control functions to check that an object has been properly
        % initialized
        flag      = hasValidInitialGuess(this)
        flag      = hasValidPrior(this)
    end
end

% ML thought 1: should perhaps the convertToOutStruct method also do some
% unit conversion, for example rescale the length unit to ,e.g., nm? 
% In that case, a scaling factor has to be supplied somehow. I would
% suggest as a construction argument, since the user presumably knows what
% units fitting will be done in, and what units one wants the final answer to be

% ML though 2: I did not set up access functions (e.g.
% this.getInitialGuess(), etc), but have no strong feelings about it.

% to use gradients: one will need to compute derivatives of the likelihood
% (diifficult for the lookup table if they are to be exact) in order to
% make use of gradients in the optimization. It would speed things up quite
% a bit I think.