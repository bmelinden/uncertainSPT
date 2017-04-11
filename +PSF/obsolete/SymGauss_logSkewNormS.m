classdef SymGauss_logSkewNormS < PSF.SymGauss
    % a symmetric Gaussian PSF model (SymGauss) where 
    %
    % 1) no background prior
    %
    % 2) the PSF width parameter lnS=log(S) has a skew-Normal distribution
    % with location parameter lnS0, std parameter lnSstd, and shape
    % parameter lnSa.    
    % see https://en.wikipedia.org/wiki/Skew_normal_distribution
    %
    % Construction:
    % P = SymGauss_logNormB('initialGuess',[mux muy lnB lnN lnS],...
    %                   'priorParameters',[lnS0 lnSstd lnSa]),
    % (the initial guess is passed on to the SymGauss constructor)
    properties (Constant)
        priorName='logSkewN psf width';
    end
    properties
        priorParameters	=[];
    end
    methods
        % Constructor
        function this = SymGauss_logSkewNormS(varargin)
            % call superclass constructor
            this@PSF.SymGauss(varargin{:});
            % sanity check on priorParameters
            if(numel(this.priorParameters)~=3)
                error('PSF.SymGauss_logSkewNormS needs priorParameters=[lnS0 lnSstd lnSa]')
            end
        end
        % compute log prior density
        function [y,dy] = logPrior(this,param)
            % y = log( p0(param) ), log (prior density)
            % dy(j) = dy / d(param(j)), gradient of y wrt param
            y=0;
            dy=zeros(size(param));
            
            % PSF width
            lnS   =param(5);
            lnS0  =this.priorParameters(1);
            lnSstd=this.priorParameters(2);
            lnSa  =this.priorParameters(3);
            
            [lnL,dlnL_dlnS]=PSF.skewGauss_logPdf(lnS,lnS0,lnSstd,lnSa);
            y=y+lnL;
            dy(5)=dlnL_dlnS;
            
        end        
    end
end