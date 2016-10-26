function [xMean,ycxMean,n]=condMeanHist(x,y,xEdges,nMin)
% [xMean,ycxMean,n]=condMeanHist(x,y,xEdges,nMin)
%
% compute conditional mean values <y|x> for x in bins specified by
% edges xEdges. Values outside the edges are ignored, and bins with too few
% counts (<nMin) give NaN.


