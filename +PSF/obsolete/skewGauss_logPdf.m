function [lnp,dlnp_dx,mx,stdx]=skewGauss_logPdf(x,mu,w,a)
% [lnp,dlnp_dx,mx,stdx]=skewGauss_logPdf(x,mu,w,a)
%
% log(pdf) for a skewed Gaussian distribution, using the definition on
% https://en.wikipedia.org/wiki/Skew_normal_distribution 
%
% x     : variable
% mu    : mean value paramerter (=mean for alpha=0).
% w     : std parameter (standard deviation when alpha=0)
% a     : skewness parameter, a=0 : standard N(mu,w)
%                             a>0 skewed towards x>mu
%                             a<0 skewed towards x<mu
%
% lnp    : log pdf
% dlnp_dx: derivative of lnp wrt x
% mx     : mean value (independent of x)
% stdx   : standard deviation (independent of x)
%
% Note that the tabulated PDF in that article is (was?) incorrect, the
% correct PDF from the main text is
%
% p(x) = 1/w/sqrt(2*pi)*exp(-0.5*(x-mu).^2/w^2)*(1+erf(a/w/sqrt(2)*x-mu));
%
% ML 20160829


y=(x-mu)/w;
lnCDFpart=log1p(erf(a/sqrt(2)*y));
lnPDFpart=-log(w)-0.5*log(2*pi)-0.5*y.^2;
lnp= lnPDFpart+lnCDFpart;
if(nargout>1)
   dlnp_dx = -(x-mu)/w^2+...
       a/w*sqrt(2/pi).*exp(-a^2*0.5*y.^2)./(1+erf(a/sqrt(2)*y));    
   % derivative verified numerically 2016-08-29, ML
end
if(nargout>2)
    delt=a*sqrt(1+a^2);
    mx=mu+w*delt*sqrt(2/pi);
    stdx=w*sqrt(1-2/pi*delt^2);
end



