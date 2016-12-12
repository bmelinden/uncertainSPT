% dat=Dest.preprocess_mixed_columns(cellData,xInd,dxInd,Lmin,centerData)
%
% assemble multi-trajectory single particle diffusion data with
% localization errors in a form for fast likelihood evaluation with
% Dest.logLlambda
%
% cellData  : a cell vector of trajectories with positions and estimated
%             errors (position posterior STANDARD DEVIATIONS). It is
%             assumed that all errors are statistically independent, i.e.,
%             correlations between errors in x(t), y(t) are neglected.
% xInd,dxInd: column indices such that
%             cellData{k}(t, xInd) = [x1 x2 ...xd]
%             cellData{k}(t,dxInd) = [e1 e2 ...ed]
%             xi : positions, ei: errors
%             Default: xInd=1, dxInd=2;
%             dxInd=[] means that the error field dat.v is left empty.
% Lmin        Minimum minimum trj length (trajectories shorter than Lmin
%             are disregarded). Optional (default 2).
% centerData: Indicator variable for position centering (optional).
%             true:  remove average position from each trajectory (default).
%             false: leave trajectory as is.
% output: data dstructure dat, with fields
% dat.x, dat.v     : positions and posterior errors (variances)
% dat.dim   : dimension
% dat.T     : vector or trajectory lengths
% dat.one, dat.end : indices of starts and ends of all trajectories.
% dat.Yone, dat.Yend: indices to start- and end-pointsof the unknown true
% trajectories (which contain one extra position in the model this data is
% constructed for). 
%
% ML 2015-03-19
%
% Note: if you want to look under the hood, note that these functions may
% use partly different index conventions than those in +EMhmm, and expect 
% precision as STANDARD DEVIATIONS, while EMhmm work with variances.


function dat=preprocess_mixed_columns(cellData,xInd,dxInd,Lmin,centerData)
%% parse input

X=cellData;
if(~exist('xInd','var'))
	xInd=1;
end
if(~exist('dxInd','var'))
	dxInd=2;
end
if(length(xInd)~=length(dxInd) && ~isempty(dxInd))
    error('ML1_preprocess_mixed_columns error: inconsistency, unequal number of position and error columns specified')
end
if(length(union(xInd,dxInd))~=length(xInd)+length(dxInd))
   warning('ML1_preprocess_mixed_columns error: possible inconsistency, some columns specified more than once.')
end
if(~exist('Lmin','var') || isempty(Lmin) )
	Lmin=2;
end
if(~exist('centerData','var') || isempty(centerData))
	centerData=true;
end
if(isempty(dxInd))
    doVariances=false;
else
    doVariances=true;
end    

%% assemble output structure
dat=struct;
dat.dim=length(xInd);
dim=dat.dim;
dat.xInd=xInd;
dat.dxInd=dxInd;

% count output size
T=zeros(size(X));
for k=1:length(X)
    T(k)=size(X{k},1);
end
% remove trajectories shorter than Lmin
X=X(T>=Lmin);
T=T(T>=Lmin);

if(isempty(T))
   error(['ML1_preprocess_mixed_columns : no trjaectories left of length >= ' int2str(Lmin)])
end
% data stacking: pack a zero-row between every trajectory to match sizes of
% data x(t) and diffusive path y(t).
dat.T=T;
dat.x=zeros(sum(T),dim);
if(doVariances)
    dat.v=zeros(sum(T),dim);
else
    dat.v=[];
end
dat.one=zeros(1,length(X),'double');
dat.end =zeros(1,length(X),'double');

ind=1;
for k=1:length(X)
    x =X{k}(:,xInd);
    if(centerData)
       x=detrend(x,'constant');
    end
    Tx=size(x,1);
    dat.one(k)=ind;    
    dat.end(k)  =ind+Tx-1;
    ind=ind+Tx;
    dat.x(dat.one(k):dat.end(k),1:dim)=x;

    if(doVariances)
        dx=X{k}(:,dxInd);
        dat.v(dat.one(k):dat.end(k),1:dim)=dx.^2;
    end
end

% indices to the hidden trajectory, which has one extra position...
% indices to hidden trajectory start- and end-points
dat.Yone=dat.one+(0:length(dat.one)-1);
dat.Yend=dat.end+(1:length(dat.one)  );



