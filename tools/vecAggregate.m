function A = vecAggregate(R,s,ind,vars)
% A = vecAggregate(R,s,ind,vars)
% aggregate cell vector or array fields from R according to 
% 1) selection indices s, or, 
% 2) if input ind is given, use ordinal indices ind, or
% 3) if neither s nor ind is given, use all elements of R.
%
% vars: cell vector of field names to aggregate. If not given, all field
% names in R are used.
%
% EAch field in every element of R is concatenated vertically in A, i.e.,
% A.foo=[R(ind(1)).foo;R(ind(2)).foo; ... ];

if(exist('ind','var') && ~isempty(ind))
    % nothing to do!
elseif(exist('s','var') && ~isempty(s))
    ind=find(s);
else
    ind=1:numel(R);
end

if(~exist('vars','var') || isempty(vars))
    vars=fieldnames(R(1));
end
if(~iscell(vars))
    vars={vars};
end
if(isempty(ind))
    error('Selection does not match any data.')
end

A=struct;
for v=1:numel(vars)
    V=R(ind(1)).(vars{v});
    for k=2:numel(ind)
        if(iscell(V) && size(V{1},1)>1)
            for m=1:numel(V)
                V{m}=[V{m};R(ind(k)).(vars{v}){m}];
            end
        elseif(size(V,1)>1)
            V=[V;R(ind(k)).(vars{v})];
        end
    end
    A.(vars{v})=V;
    %disp(vars{v})
end
