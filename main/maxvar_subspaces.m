%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subspaces
% Copyright (C) 2016 Gamaleldin F. Elsayed and John P. Cunningham
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gamaleldin F. Elsayed
% 
% maxvar_subspaces.m
%
%     Inputs: for the ith subspace
%           DataStruct(i).A : contains sample data  (#condition . #times X 
%                               #neurons) to identify its subspace;          
%           DataStruct(i).dim : dimensionality of this subspace;
%     Optional Inputs:
%           DataStruct(i).bias : weight of the ith data relative to the
%           others.
%     Outputs: 
%           QSubspaces(i).Q : orthonormal vectors of the identified
%           ith subspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples
% data1 = randn(400,2)*orth(randn(100, 2))'+0.1*randn(400,100); % gen data
% data2 = randn(400,4)*orth(randn(100, 4))'+0.1*randn(400,100); % gen data
% DataStruct(1).A = data1;
% DataStruct(1).dim = 2;
% DataStruct(2).A = data2;
% DataStruct(2).dim = 4;
% [QSubspaces] = maxvar_subspaces(DataStruct)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [QSubspaces] = maxvar_subspaces(DataStruct)
figFlg = false;
numSubspaces = length(DataStruct);
totalDim = 0;
bias = ones(numSubspaces,1);
if isfield(DataStruct, 'bias')
bias = [DataStruct.bias];
end
bias = bias./sum(bias);
Q = [];

for j = 1:numSubspaces
    Cj = cov(DataStruct(j).A); 
    d = DataStruct(j).dim; 
    [V, S] = svd(Cj);
    S = diag(S);
    covM(:,:,j) = Cj; 
    dim(j) = d;
    totalDim = totalDim+d;
%     Q = [Q V(:,end-d+1:end)];
    Q = [Q V(:,1:d)];
    normFact(j) = bias(j)./sum(S(1:d));
end

[Q1,~, Q2] = svd(Q, 'econ');
Q = Q1*Q2';

maxIter = 10000;
tic;
[ Q , ~ , info] = minimize_stiefel_trust( 'maxvar_subspaces_objfn' , Q , maxIter , covM, dim, normFact );
iterTime = toc;
f = 1-[info.cost];




[~, ~, QSubspaces] = maxvar_subspaces_objfn(Q, covM, dim, normFact);


%% Align the axes within a subspace by variance high to low
for j = 1:numSubspaces
    Aj = DataStruct(j).A;
    Qj = QSubspaces(j).Q;
    Aj_proj = bsxfun(@minus, Aj, mean(Aj))*Qj;
    [V, ~, ~] = svd(cov(Aj_proj));
    Qj = Qj*V;
    QSubspaces(j).Q = Qj; % rank top to low variance    
end

%%
fprintf('Total iterations: %5d of %d,      Final Obj. Value: %g,  Total time: (%.5f sec)', length(f),  maxIter, f(end), sum(iterTime));
fprintf('\n');


if figFlg
    figure; 
    hold on
    plot(ones(length(f),1),'k:')
    plot(1:length(f), f,'b.-')
    xlabel('iter')
    ylabel('Separation index')
    set(gca, 'FontSize', 16)
    ylim([0 1.1])
    set(gca,'Ytick', (0:0.2:1))
    hold off
end
end


