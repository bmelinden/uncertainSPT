classdef SymGaussS0 < PSF.PSFmodel
    % symmetric Gaussian PSF model with width offset, described by
    % parameters 
    % p = [muX muY lnB lnN lndS], and with density
    %
    % E = B + N/2/pi/S^2*exp(-0.5*((xx-muX)/S)^2-0.5*((yy-muY)/S)^2), 
    %
    % where B= exp(lnB), N = exp(lnN), S = S0*(1+exp(lndS)), 
    % and the PSF width scale is S0=0.21*lambda/NA. 
    % The offset is the Gaussian approximation of a perfectly focused spot,
    % and therefore a physics-based lower on the PSF width. Note that while
    % lambda and NA are object properties, only the S0 property is used
    % in actual computations.
    %
    % Note that this is a continuous model, so that exp(lnB) is the
    % background intensity per unit area, not per pixel. (But when doing
    % math in pixel units the area per pixel is 1, and then this does not
    % matter.)

    properties (Constant)
        modelName = 'Symmetric Gaussian s_min';
    end
    properties
       initialGuess=[]; 
        lambda=[];
        NA=[];  
        S0=[];
    end
    
    methods (Access = public)
        % Constructor
        function this = SymGaussS0(varargin)            
            this@PSF.PSFmodel(varargin{:});
            % sanity check for initialGuess
            if(numel(this.initialGuess)~=5)
                error('SymGaussS0 needs 5 initialGuess elements')
            end
            if( ~isempty(this.NA) && ~isempty(this.lambda))
                this.S0=0.21*this.lambda/this.NA;
            elseif(isempty(this.S0))
                error('SymGaussS0 must be initialized with either NA,lambda or S0.')
            end
        end
        % PSFmodel functions in this class level
        function [E,dEdp]= psfDensity(this, xx, yy, param)
        %   [E,dEdp]= psfDensity(this, xx, yy, param)  
            muX=param(1);
            muY=param(2);
            B=exp(param(3));
            N=exp(param(4));
            S=this.S0*(1+exp(param(5)));
            
            S2=S^2;            
            NEexp=N/S2/2/pi*exp(-1/2/S2*((muX-xx).^2+(muY-yy).^2));
            E=B+NEexp;

            if(nargout>1)
                dEdp=zeros(size(E,1),size(E,2),5);
                dEdp(:,:,1)=-(muX-xx).*NEexp/S2; % dE_dmuX 
                dEdp(:,:,2)=-(muY-yy).*NEexp/S2; % dE_dmuY 
                dEdp(:,:,3)=ones(size(E))*B; % dE_dlnBG
                dEdp(:,:,4)= NEexp; % dE_dlnN 
                %dEdp(:,:,5)= NEexp.*(-2+((muX-xx).^2+(muY-yy).^2)/S2); % dE_dlndS 
                error('derivative wrt lndS not implemented yet')
                % compute parameter dependent pixel intensities for symmetric Gaussian, with derivatives
                % ML 2015-11-17 : partial derivatives validated, except for
                % the new offset width.
            end
        end
        function [outStruct,outPar] = convertToOutStruct(this,inPar, inStruct)
            % [outStruct,outPar] = convertToOutStruct(this,inPar, inStruct)
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
            % outStruct.std         : spot width, stdandard deviation   = S0*(1+exp(lndS))
            % outPar                : [B N S]
            
            if(nargin==3)
                outStruct = inStruct;
            else
                outStruct=struct;
            end
            
            B=exp(inPar(3));
            N=exp(inPar(4));
            S=this.S0*(1+exp(inPar(5)));

            outStruct.background=B;
            outStruct.amplitude =N;
            outStruct.std       =S;
            outPar=[B N S];
        end        
    end
end
