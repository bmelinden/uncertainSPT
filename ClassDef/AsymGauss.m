classdef AsymGauss < PSFModel
    
    properties (Constant)
        % asymmetric Gaussian PSF model and partial derivatives, 
        % parameterized by parameters
        % p=[ muX muY lnB lnN lnS1 lnS2 v] :
        % E = exp(lnB) + f, with 
        % f = N/2/pi/S1/S2*exp(-0.5*( cos(v)*(xx-muX)+sin(v)*(yy-muY))^2/S1^2
        %                      -0.5*(-sin(v)*(xx-muX)+cos(v)*(yy-muY))^2/S2^2)
        %
        % Note that this is a continuous model, so that exp(lnB) is the background
        % intensity per unit area, not per pixel. But when doing math in pixel
        % units the area per pixel is 1, and then this does not matter.
        modelName = 'Asymmetric Gaussian';
        
    end
    
    methods
        % Constructor
        function this = AsymGauss(varargin)
            if rem(length(varargin),2) == 1
                varargin = varargin{1};
            end
            for i = 1 : 2 : length(varargin)
                this.(varargin{i}) = varargin{i+1};
            end
            
            if isempty(this.parameters) && ~isempty(this.initialGuess)
                this.parameters = this.initialGuess;
            end
        end
        
        function [E,dE_dmuX,dE_dmuY,dE_dlnB,dE_dlnN,dE_dlnS1,dE_dlnS2,dE_dv] = runModel(this, xx, yy, param)
            
            [muX, muY, lnB, lnN, lnS1, lnS2, v] = translateFitParameters(param);

            sig1=exp(lnS1);
            sig2=exp(lnS2);
            c =cos(v);
            s =sin(v);
            dx=xx-muX;
            dy=yy-muY;

            NEexp=1/2/pi*exp(lnN-lnS1-lnS2-1/2*(((c*dx+s*dy)/sig1).^2+((-s*dx+c*dy)/sig2).^2));
            E=NEexp+exp(lnB);

            if(nargout>1)
                dE_dmuX =NEexp.*( c/sig1^2*(c*dx+s*dy)-s/sig2^2*(-s*dx+c*dy));
                dE_dmuY =NEexp.*( s/sig1^2*(c*dx+s*dy)+c/sig2^2*(-s*dx+c*dy));
                dE_dlnB = ones(size(E))*exp(lnB);
                dE_dlnN = NEexp;
                dE_dlnS1 = NEexp.*(-1+((c*dx+s*dy)/sig1)^2);
                dE_dlnS2 = NEexp.*(-1+((-s*dx+c*dy)/sig2)^2);
                dE_dv = NEexp.*(c*dx+s*dy).*(s*dx-c*dy)*(1/sig1-1/sig2)*(1/sig1+1/sig2);
                % compute parameter dependent pixel intensities for symmetric Gaussian, with derivatives
            end
        end
        
        function [outStruct] = convertToOutStruct(this, inStruct, inPar)
            % Convert psf fit parameters to parameter struct, using a symmetric
            % Gaussian PSF model, parameterized by 
            % p=[ muX muY lnB lnN lnS] :
            % E = exp(lnB) + f, with 
            % f = N/2/pi/S^2*exp(-0.5*((xx-muX)/S)^2-0.5*((yy-muY)/S)^2)
            %
            % Note that this is a continuous model, so that exp(lnB) is the background
            % intensity per unit area, not per pixel. But when doing math in pixel
            % units the area per pixel is 1, and then this does not matter.
            %
            % If no inPar input, the currently stored values for the
            % parameters are returned.
            %
            % output: parameter struct outStruct, with fields
            % outStruct.background  : background intensity photons/area = exp(lnB)
            % outStruct.amplitude   : spot amplitude (photons)          = exp(lnN)
            % outStruct.std         : spot width, stdandard deviation   = exp(lnS)
            
            outStruct = inStruct;
            
            if nargin > 2
                outStruct.background=exp(inPar(3));
                outStruct.amplitude =exp(inPar(4));
                outStruct.std1      =exp(inPar(5));
                outStruct.std2      =exp(inPar(6));
                outStruct.angle     =inPar(7);
            else
                outStruct.background=exp(this.parameters(3));
                outStruct.amplitude =exp(this.parameters(4));
                outStruct.std1      =exp(this.parameters(5));
                outStruct.std2      =exp(this.parameters(6));
                outStruct.angle     =this.parameters(7);
            end
        end
        
        function [muX, muY, lnB, lnN, lnS1, lnS2, v] = translateFitParameters(this, param)
            muX=param(1);
            muY=param(2);
            lnB=param(3);
            lnN=param(4);
            lnS1=param(5);
            lnS2=param(6);
            v=param(7);
        end
    end
end