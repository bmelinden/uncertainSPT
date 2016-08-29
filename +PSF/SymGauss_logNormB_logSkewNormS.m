classdef SymGauss_logNormB_logSkewNormS < PSF.SymGauss
    % a symmetric Gaussian PSF model (SymGauss) where 
    %
    % 1) the background intensity B has a log-normal prior with location
    % parameter lnB0 and scale parameter lnBstd, which means that
    % lnB is N(lnB0,lnBstd^2), or 
    % <lnB> = lnB0, Std(lnB) =lnBstd, <B> = <exp(lnB)> = exp(lnB0+lnBstd^2/2),
    % std(B) = exp(lnB0+lnBstd^2/2)*sqrt(exp(lnBstd^2)-1),
    % see https://en.wikipedia.org/wiki/Log-normal_distribution.
    %
    % 2) the PSF width parameter lnS=log(S) has a skew-Normal distribution
    % with location parameter lnS0, std parameter lnSstd, and shape
    % parameter lnSa.    
    % see https://en.wikipedia.org/wiki/Skew_normal_distribution
    %
    % Construction:
    % P = SymGauss_logNormB('initialGuess',[mux muy lnB lnN lnS],...
    %                   'priorParameters',[lnB0 lnBstd lnS0 lnSstd lnSa]),
    % (the initial guess is passed on to the SymGauss constructor)
    properties (Constant)
        priorName='logN background, logSkewN psf width';
    end
    properties
        priorParameters	=[];
    end
    methods
        % Constructor
        function this = SymGauss_logNormB_logSkewNormS(varargin)
            % call superclass constructor
            this@PSF.SymGauss(varargin{:});
            % sanity check on priorParameters
            if(numel(this.priorParameters)~=5)
                error('PSF.SymGauss_logNormB_logSkewNormS needs priorParameters=[lnB0 lnBstd lnS0 lnSstd lnSa]')
            end
        end
        % compute log prior density
        function [y,dy] = logPrior(this,param)
            % y = log( p0(param) ), log (prior density)
            % dy(j) = dy / d(param(j)), gradient of y wrt param
            y=0;
            dy=zeros(size(param));
            
            % background 
            lnB=param(3);
            lnB0=this.priorParameters(1);
            lnBstd=this.priorParameters(2);
           
            y = y - 1/2 *(lnB - lnB0).^2./(lnBstd^2);
            dy(3)=- (lnB - lnB0)./(lnBstd^2);
           
            % PSF width
            lnS   =param(5);
            lnS0  =this.priorParameters(3);
            lnSstd=this.priorParameters(4);
            lnSa  =this.priorParameters(5);
            
            [lnL,dlnL_dlnS]=skewGauss_logPdf(lnS,lnS0,lnSstd,lnSa);
            y=y+lnL;
            dy(5)=dlnL_dlnS;
            
        end        
    end
end