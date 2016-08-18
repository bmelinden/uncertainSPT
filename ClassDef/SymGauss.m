classdef SymGauss < PSFModel
    
    properties (Constant)
        % symmetric Gaussian PSF model and partial derivatives, parameterized by
        % parameters = [muX muY lnB lnN lnS] :
        % E = exp(lnB) + f, with 
        % f = N/2/pi/S^2*exp(-0.5*((xx-muX)/S)^2-0.5*((yy-muY)/S)^2)
        modelName = 'Symmetric Gaussian';
    end
    
    methods
        % Constructor
        function this = SymGauss(varargin)

        end
        
        function [E,dE_dmuX,dE_dmuY,dE_dlnBG,dE_dlnN,dE_dlnS] = runModel(this, xx, yy, param)
            
            [muX, muY, lnB, N, S2] = translateFitParameters(param);
            
            NEexp=N/S2*exp(-1/2/S2*((muX-xx).^2+(muY-yy).^2))/2/pi;

            E=exp(lnB)+NEexp;

            if(nargout>1)
                dE_dmuX =-(muX-xx).*NEexp/S2;
                dE_dmuY =-(muY-yy).*NEexp/S2;
                dE_dlnBG=ones(size(E))*exp(lnB);
                dE_dlnN = NEexp;
                dE_dlnS = NEexp.*(-2+((muX-xx).^2+(muY-yy).^2)/S2);
                % compute parameter dependent pixel intensities for symmetric Gaussian, with derivatives
                % ML 2015-11-17 : partial derivatives validated
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
            
            outStruct.background=exp(inPar(3));
            outStruct.amplitude =exp(inPar(4));
            outStruct.std       =exp(inPar(5));
        end
        
        function [muX, muY, lnB, N, S2] = translateFitParameters(this, param)
            muX=param(1);
            muY=param(2);
            lnB=param(3);
            N=exp(param(4));
            S2=exp(2*param(5)); % S^2
        end
        
    end
    
end