function [lnp,mx,stdx]=skewGauss_logPdf(x,mu,w,a)
% lnp=skewGauss_logPdf(x,mu,sigma,alpha)
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
% lnp   : log pdf
% mx    : mean value (independent of x)
% stdx  : standard deviation (independent of x)
%
% ML 20160509


y=(x-mu)/w;
lnCDFpart=log1p(erf(a*y/sqrt(2)));
lnPDFpart=-log(w)-0.5*log(2*pi)-0.5*y.^2;
lnp= lnPDFpart+lnCDFpart;

if(nargout>1)
    delt=a*sqrt(1+a^2);
    mx=mu+w*delt*sqrt(2/pi);
    stdx=w*sqrt(1-2/pi*delt^2);
end