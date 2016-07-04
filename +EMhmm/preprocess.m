function dat=preprocess(varargin)
% dat=EMhmm.preprocess(runinput)
% dat=EMhmm.preprocess(runinput,dim)
% dat=EMhmm.preprocess(X,varX,dim)
% dat=EMhmm.preprocess(X,varX,dim,misc)
%
% Assemble single particle diffusion data for diffusive HMM analysis
% runinput  : a runinput file or runinput structure. This functionality is
% not part of the EMhmm package.
% dim       : number of data dimensions (x,y,z,...)
% X         : cell vector of position trajectories, or single trajectory
% varX      : cell vector of position uncertainties (posterior variances),
%             or a vector of such variances.
% misc      : cell vector of some other field one would like to keep track
%             of. The row elements of each cell element are organized in
%             the same way as the positions, for easy comparison.
% ML 2016-05-13

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocess, data preprocessor for diffusive HMM analysis
% =========================================================================
% 
% Copyright (C) 2016 Martin Lindén
% 
% E-mail: bmelinden@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by linking or combining it
%  with Matlab or any Matlab toolbox, the licensors of this Program grant you 
%  additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% parse input
% if an existing file, generate options structure
doMisc=false;
if(nargin<3)
    runinput=varargin{1};
    if(ischar(runinput) && exist(runinput, 'file')==2)
        runinputfile = runinput;
        opt=VB3_getOptions(runinputfile);
        X=VB6_readData(opt);
    elseif(isstruct(runinput)) % if already an options structure
        opt=runinput;
        %runinputfile=opt.runinputfile;
        X=VB6_readData(opt);
    end
    if(nargin==2)
        dim=varargin{2};
    else
        dim=opt.dim;
    end
    clear opt runinput runinputfile    
elseif(nargin>=3)
    X=varargin{1};
    varX=varargin{2};
    dim=varargin{3};    
    if(nargin==4)
        doMisc=true;
        misc=varargin{4};
    end
else
    error(['Not a valid input, aborting SPT_preprocess']);
end
if(~iscell(X) && ~iscell(varX))
    X={X};
    varX={varX};
    if(doMisc && ~iscell(misc))
        misc={misc};
    end
end
%% assemble output structure
dat=struct;
dat.dim=dim;

% count output size
T=zeros(size(X));
for k=1:length(X)
    T(k)=size(X{k},1);
end
if(~isempty(find(T<2,1)))
   warning('EMhmm.preprocess: data contains traces with no steps. FIX!')
end
% data stacking: pack a zero-row between every trajectory to match sizes of
% data x(t) and diffusive path y(t).
dat.T=T;
dat.x=zeros(sum(T+1),dim);
dat.v=zeros(sum(T+1),dim);
dat.one=zeros(1,length(X),'double');
dat.end =zeros(1,length(X),'double');
if(doMisc)
    miscColumns=size(misc{1},2);
   dat.misc=zeros(sum(T+1),miscColumns);
end


ind=1;
for k=1:length(X)
    x=X{k}(:,1:dim);
    v=varX{k}(:,1:dim);
    Tx=size(x,1);
    % a sanity check: v>0
    if(~isempty(find(isfinite(v(:)).*(v(:)<=0))))
       error(['EMhmm.preprocess found variance<=0 in trj ' int2str(k) ])
    end    
    dat.one(k)=ind;
    dat.end(k)=ind+Tx-1;
    ind=ind+Tx;
    dat.x(dat.one(k):dat.end(k),1:dim)=x; 
    dat.v(dat.one(k):dat.end(k),1:dim)=v;
    if(doMisc)
       dat.misc( dat.one(k):dat.end(k),1:miscColumns)=misc{k};
    end
    ind=ind+1;
end
dat.x(isnan(dat.x))=0;

