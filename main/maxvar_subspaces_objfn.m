%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subspaces
% Copyright (C) 2016 Gamaleldin F. Elsayed and John P. Cunningham
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gamaleldin F. Elsayed

function [f, gradQ, QSubspaces] = maxvar_subspaces_objfn(Q, C, dim, normFact)
numSubspaces = length(dim);
gradQ = nan(size(Q));
dim = [0;dim(:)];
f = nan(numSubspaces,1);
QSubspaces = struct([]);
for j = 1:numSubspaces
    Qj = Q(:,sum(dim(1:j))+1:sum(dim(1:j))+dim(j+1));
    Cj = C(:,:,j);
    normFactj = normFact(j);
    gradQ(:,sum(dim(1:j))+1:sum(dim(1:j))+dim(j+1)) = Cj*Qj*normFactj;
    f(j) = normFactj*trace(Qj'*Cj*Qj);
    QSubspaces(j).Q = Qj;
    QSubspaces(j).costFn = f(j);
end
f = 1-sum(f);
gradQ = -gradQ;
end