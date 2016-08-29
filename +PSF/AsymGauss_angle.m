classdef AsymGauss_angle < PSF.PSFModel
    % Asymmetric Gaussian PSF model described by parameters
    % p=[ muX muY lnB lnN lnS1 lnS2 v], and with density
    %
    % E = B + N/2/pi/S1/S2*exp(
    %           -0.5*( cos(v)*(xx-muX)+sin(v)*(yy-muY))^2/S1^2
    %           -0.5*(-sin(v)*(xx-muX)+cos(v)*(yy-muY))^2/S2^2),
    %
    % where B=exp(lnB), N=exp(lnN), S1=exp(lnS1), S2=exp(lnS2), and v
    % is an angle in radians. The model can also compute derivatives wrt
    % the parameters.
    %
    % Note that this is a continuous model, so that exp(lnB) is the
    % background intensity per unit area, not per pixel. (However, when
    % doing math in pixel units the area per pixel is 1, and then this does
    % not matter.)
    
    properties (Constant)
        modelName = 'Asymmetric Gaussian, std12,angle parameterization';
    end
    properties
        initialGuess=[];
    end
    methods (Access = public)
        % Constructor, requires a name-value pair
        % 'initialGuess',[ muX muY lnB lnN lnS1 lnS2 v ]
        function this= AsymGauss_angle(varargin)
            this@PSF.PSFModel(varargin{:});
            % sanity check for initialGuess
            if(numel(this.initialGuess)~=7)
                error('AsymGauss_angle needs 7 initialGuess elements')
            end
        end
        
        % PSFmodel functions in this class level
        function [E,dEdp]= psfDensity(this, xx, yy, param)
            
            muX=param(1);
            muY=param(2);
            lnB=param(3);
            lnN=param(4);
            lnS1=param(5);
            lnS2=param(6);
            v=param(7);
            
            sig1=exp(lnS1);
            sig2=exp(lnS2);
            c =cos(v);
            s =sin(v);
            dx=xx-muX;
            dy=yy-muY;
            
            NEexp=1/2/pi*exp(lnN-lnS1-lnS2-1/2*(((c*dx+s*dy)/sig1).^2+((-s*dx+c*dy)/sig2).^2));
            E=NEexp+exp(lnB);
            
            % ML: derivatives verified numerically 2016-05-31
            if(nargout>1)
                dEdp=zeros(size(E,1),size(E,2),7);
                dEdp(:,:,1)=NEexp.*( c/sig1^2*(c*dx+s*dy)-s/sig2^2*(-s*dx+c*dy));% dE_dmuX
                dEdp(:,:,2)=NEexp.*( s/sig1^2*(c*dx+s*dy)+c/sig2^2*(-s*dx+c*dy));% dE_dmuY
                dEdp(:,:,3)= ones(size(E))*exp(lnB); % dE_dlnB
                dEdp(:,:,4)= NEexp; % dE_dlnN
                dEdp(:,:,5)= NEexp.*(-1+((c*dx+s*dy)/sig1)^2); % dE_dlnS1
                dEdp(:,:,6)= NEexp.*(-1+((-s*dx+c*dy)/sig2)^2);% dE_dlnS2
                dEdp(:,:,7)= NEexp.*(c*dx+s*dy).*(s*dx-c*dy)*(1/sig1-1/sig2)*(1/sig1+1/sig2);% dE_dv
                % compute parameter dependent pixel intensities for symmetric Gaussian, with derivatives
            end
        end
        
        % convert a parameter set to a nice-looking parameter struct
        function outStruct = convertToOutStruct(this,inPar, inStruct)
            % outStruct = convertToOutStruct(inPar, inStruct)
            %
            % Convert psf fit parameters to a nice-looking parameter
            % struct
            %
            % inPar     : parameter vector to convert.
            % inStruct  : struct to which parameter fields are (give an
            %             empty struct if no struct is at hand).
            %
            % output: parameter struct outStruct, with fields
            % outStruct.background  : background intensity photons/area = exp(lnB)
            % outStruct.amplitude   : spot amplitude (photons)          = exp(lnN)
            % outStruct.std12       : spot principal widths             = exp([lnS lnS2])
            % outStruct.angle       : spot orientation                  = v
            
            outStruct = inStruct;
            
            B=exp(inPar(3));
            N=exp(inPar(4));
            S1=exp(inPar(5));
            S2=exp(inPar(6));
            v=inPar(7);
            
            % make the angle refer to the largest principal direction?
            %%% ML: good idea, but leave for later
            
            outStruct.background=B;
            outStruct.amplitude =N;
            outStruct.std12     =[S1 S2];
            outStruct.angle=v;
        end
    end
end
