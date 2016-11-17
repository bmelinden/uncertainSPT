function [FcMean,xcMean,n,dx]=condMeanHist(x,y,xE,F,nMin,setNaN)
% [FcMean,xcMean,n,dx]=condMeanHist(x,y,xE,F,nMin,setNaN)
%
% compute conditional mean values <F(y)|x in bin j> for x in bins specified by
% edges xE. Values outside the edges are ignored, and bins with too few
% counts (<nMin) give NaN.
%
%
% x     : independent data (for binning)
% y     : dependent data (for averaging)
% xE    : bin edges for x, entries with x-values outside these edges are
%         ignored
% F     : function to average (function handle or anonymous function).
%         Default: @std
% nMin  : optional lower bound on number of counts per bin. 
%         Default: 1.
% setNaN: if true, bins with n<nMin give ycMean=NaN. If false, such bins
%         are removed from the output. 
%         Default: true.
%
% FcMean : conditional function values, Fc(j) = F(y | xE(j)<= x < xE(j+1)), 
%         j=1,2,..., length(xE)-1. Note that it must retur a single value
%         from multiple inputs. Default: @mean
% xcMean: conditional mean value xcMean(j) = <x| xE(j)<= x < xE(j+1)>. If
%         n<nMin, the arithmetic mean (xE(j)+xE(j+1))/2 is returned.
% n     : bin counts 
% dx    : bin width, dx(j) =  xE(j+1)-xE(j).
% 
% ML 2016-10-26

% parameter check
if(~exist('F','var') || isempty(F))
    F=@mean;
end

if(~exist('nMin','var') || isempty(nMin))
    nMin=1;
end
if(~exist('setNaN','var') || isempty(setNaN))
    setNaN=true;
end

dx=diff(xE);
if( ~isempty(find(dx<=0,1)))
   error('condMeanHist requires edges xE to be strictly increasing') 
end
% construct histograms and conditional 
FcMean = zeros(size(dx));
xcMean = zeros(size(dx));
n      = zeros(size(dx));
NB=length(n);

for k=1:NB
    ind = xE(k) <= x & x < xE(k+1);
    nk=sum(ind);
    if(nk>=nMin)
        xcMean(k)=mean(x(ind));
        FcMean(k)=F(y(ind));
        n(k)     =nk;
    end
end

ind= n<nMin;
binMid=(xE(1:end-1)+xE(2:end))/2;
xcMean(ind)=binMid(ind);
if(setNaN)
    FcMean(ind)=nan;
else
    xcMean=xcMean(~ind);
    FcMean=FcMean(~ind);
    n = n(~ind);
    dx=dx(~ind);
end
