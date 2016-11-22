classdef AsymGaussS0 < PSF.PSFmodel
    % Asymmetric Gaussian PSF model with PSF width offset, described by
    % parameters p=[ muX muY lnB lnN lndS1 lndS2 v], and spot density
    %
    % E = B + N/2/pi/S1/S2*exp(
    %           -0.5*( cos(v)*(xx-muX)+sin(v)*(yy-muY))^2/S1^2
    %           -0.5*(-sin(v)*(xx-muX)+cos(v)*(yy-muY))^2/S2^2),
    %
    % where B=exp(lnB), N=exp(lnN), S1=S0*(1+exp(lnS1)),
    % S2=S0*(1+exp(lnS2)), and v is an angle in radians. The minimum PSF
    % width S0 is either given directly, or in terms of lambda and NA,
    % S0=0.21*lambda/NA. This parameterization extends that of SymGaussS0
    % to asymmetric spots.
    %
    % The model can also compute derivatives wrt the parameters.
    %
    % Note that this is a continuous model, so that exp(lnB) is the
    % background intensity per unit area, not per pixel. (However, when 
    % doing math in pixel units the area per pixel is 1, and then this does
    % not matter.)
    
    properties (Constant)
        modelName = 'Asymmetric GaussianS0';
    end
    properties
        initialGuess=[];
        lambda=[];
        NA=[];
        S0=[];
    end
    methods (Access = public)
        % Constructor, requires a name-value pair
        % 'initialGuess',[ muX muY lnB lnN lnS1 lnS2 v ]
        function this= AsymGaussS0(varargin)
            this@PSF.PSFmodel(varargin{:});
            % sanity check for initialGuess
            if(numel(this.initialGuess)~=7)
                error('AsymGaussS0 needs 7 initialGuess elements')
            end
            if( ~isempty(this.NA) && ~isempty(this.lambda))
                this.S0=0.21*this.lambda/this.NA;
            elseif(isempty(this.S0))
                error('AsymGaussS0 must be initialized with either NA,lambda or S0.')
            end

        end
        
        % PSFmodel functions in this class level
        function [E,dEdp]= psfDensity(this, xx, yy, param)
            
            muX=param(1);
            muY=param(2);
            lnB=param(3);
            lnN=param(4);
            lndS1=param(5);
            lndS2=param(6);
            v=param(7);
            
            sig0=this.S0;
            sig1=sig0*(1+exp(lndS1));
            sig2=sig0*(1+exp(lndS2));
            c =cos(v);
            s =sin(v);
            dx=xx-muX;
            dy=yy-muY;
            
            NEexp=1/2/pi/sig1/sig2*exp(lnN-1/2*(((c*dx+s*dy)/sig1).^2+((-s*dx+c*dy)/sig2).^2));
            E=NEexp+exp(lnB);
            
            % ML: derivatives verified numerically 2016-05-31
            % ML: derivatives wrt lndS1,2 verified numerically 2016-11-21
            if(nargout>1)
                dEdp=zeros(size(E,1),size(E,2),7);
                dEdp(:,:,1)=NEexp.*( c/sig1^2*(c*dx+s*dy)-s/sig2^2*(-s*dx+c*dy));% dE_dmuX
                dEdp(:,:,2)=NEexp.*( s/sig1^2*(c*dx+s*dy)+c/sig2^2*(-s*dx+c*dy));% dE_dmuY
                dEdp(:,:,3)= ones(size(E))*exp(lnB); % dE_dlnB
                dEdp(:,:,4)= NEexp; % dE_dlnN
                dEdp(:,:,5)= NEexp.*( (sig1-sig0)/sig1*(-1+(( c*dx+s*dy)/sig1)^2)); % dE_dlnS1
                dEdp(:,:,6)= NEexp.*( (sig2-sig0)/sig2*(-1+((-s*dx+c*dy)/sig2)^2)); % dE_dlnS2
                dEdp(:,:,7)= NEexp.*(c*dx+s*dy).*(s*dx-c*dy)*(1/sig1-1/sig2)*(1/sig1+1/sig2);% dE_dv
                % compute parameter dependent pixel intensities for symmetric Gaussian, with derivatives
            end
        end
        
        % convert a parameter set to a nice-looking parameter struct
        function [outStruct,outPar] = convertToOutStruct(this,inPar, inStruct)
            % [outStruct,outPar] = convertToOutStruct(inPar, inStruct)
            %
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
            % outStruct.std1        : spot principal width 1            = exp(lnS1)
            % outStruct.std2        : spot principal width 2            = exp(lnS2)
            % outStruct.angle       : spot orientation                  = v
            % outPar                : [B N S1 S2 v]
            if(nargin==3)
                outStruct = inStruct;
            else
                outStruct=struct;
            end
            
            B=exp(inPar(3));
            N=exp(inPar(4));
            S1=this.S0*(1+exp(inPar(5)));
            S2=this.S0*(1+exp(inPar(6)));
            v=inPar(7);
            
            % make the angle refer to the largest principal direction?
            %%% ML: good idea, but leave for later
            
            outStruct.background=B;
            outStruct.amplitude =N;
            outStruct.std1      =S1;
            outStruct.std2      =S2;
            outStruct.angle=v;
            
            outPar=[B N S1 S2 v];
        end
    end
end
