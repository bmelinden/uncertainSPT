function Y=listToCells(dat,X)
% split a vector reprepresentation, such as W.Y.myU or the sMaxP or Viterbi
% state output, into a cell vector, where the cells have the same size as
% the input trajectories.
Y=cell(1,numel(dat.i0));
for k=1:numel(Y)
    Y{k}=X(dat.i0(k):dat.i1(k),:);
end

