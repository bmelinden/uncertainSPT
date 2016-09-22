classdef SymGauss < PSF.PSFmodel
    % symmetric Gaussian PSF model, described by parameters
    % p = [muX muY lnB lnN lnS], and with density
    %
    % E = B + N/2/pi/S^2*exp(-0.5*((xx-muX)/S)^2-0.5*((yy-muY)/S)^2), 
    %
    % where B= exp(lnB), N = exp(lnN), S = exp(lnS).
    %
    % Note that this is a continuous model, so that exp(lnB) is the
    % background intensity per unit area, not per pixel. (But when doing
    % math in pixel units the area per pixel is 1, and then this does not
    % matter.)

    properties (Constant)
        modelName = 'Symmetric Gaussian';
    end
    properties
       initialGuess=[]; 
    end
    
    methods (Access = public)
        % Constructor
        function this = SymGauss(varargin)
            this@PSF.PSFmodel(varargin{:});
            % sanity check for initialGuess
            if(numel(this.initialGuess)~=5)
                error('SymGauss needs 5 initialGuess elements')
            end
        end
        % PSFmodel functions in this class level
        function [E,dEdp]= psfDensity(this, xx, yy, param)
        %   [E,dEdp]= psfDensity(this, xx, yy, param)  
            [muX, muY, B, N, S] = this.translateFitParameters(param);
            S2=S^2;            
            NEexp=N/S2/2/pi*exp(-1/2/S2*((muX-xx).^2+(muY-yy).^2));
            E=B+NEexp;

            if(nargout>1)
                dEdp=zeros(size(E,1),size(E,2),5);
                dEdp(:,:,1)=-(muX-xx).*NEexp/S2; % dE_dmuX 
                dEdp(:,:,2)=-(muY-yy).*NEexp/S2; % dE_dmuY 
                dEdp(:,:,3)=ones(size(E))*B; % dE_dlnBG
                dEdp(:,:,4)= NEexp; % dE_dlnN 
                dEdp(:,:,5)= NEexp.*(-2+((muX-xx).^2+(muY-yy).^2)/S2); % dE_dlnS 
                % compute parameter dependent pixel intensities for symmetric Gaussian, with derivatives
                % ML 2015-11-17 : partial derivatives validated
            end
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
            % outStruct.amplitude   : spot amplitude (photons)          = exp(lnN)
            % outStruct.std         : spot width, stdandard deviation   = exp(lnS)
            
            if(nargin==3)
                outStruct = inStruct;
            else
                outStruct=struct;
            end
            
            [~, ~, B, N, S] = this.translateFitParameters(inPar);
            
            outStruct.background=B;
            outStruct.amplitude =N;
            outStruct.std       =S;
        end        
    end
    % protected methods
    methods (Access = protected)
        function [muX, muY, B, N, S] = translateFitParameters(~,param)
            muX=param(1);
            muY=param(2);
            B=exp(param(3));
            N=exp(param(4));
            S=exp(param(5));
        end
        
    end    
end
