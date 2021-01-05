%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <subspaces>
% Copyright (C) 2016 Gamaleldin F. Elsayed and John P. Cunningham 
%       (see full notice in README)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demonstration of how to use this code package 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

startup
% define spaces


dim = 2;
theta = [0 45 90];


for i = 1:3
A = randn(1000,dim);
Q = orth(randn(dim));
Q1 = Q(:,1);
rotMtx = [cos(theta(i)*pi/180) sin(theta(i)*pi/180); 
        - sin(theta(i)*pi/180) cos(theta(i)*pi/180)];
Q2 = rotMtx*Q(:,2);
    
A1 = A*(Q1*Q1')+0.05*randn(size(A)); %  generate data at one subspace
A2 = 0.25*(A*(Q2*Q2')+0.05*randn(size(A))); % generate data at another subspace with different variance
DataStruct(1).A = A1; % epoch 1 data 
DataStruct(1).dim = 1; % dimensionality of subspace 1
DataStruct(2).A = A2; % epoch 2 data 
DataStruct(2).dim = 1; % dimensionality of subspace 2
[QSubspaces] = maxvar_subspaces(DataStruct);
Qi1 = QSubspaces(1).Q; % identified orthonormal basis for subspace 1
Qi2 = QSubspaces(2).Q; % identified orthonormal basis for subspace 2

hf(i) = figure;hold on
set(gca,'plotboxaspectratio',[1 1 1])
set(hf(i), 'color', [1 1 1]);
box on
axis(2*[-1 1 -1 1]) 
set(gca,'xtick', [])
set(gca,'ytick', [])
hold on
plot(A1(:,1), A1(:,2), 'r.', 'markersize',4)
plot(A2(:,1), A2(:,2), 'go', 'markersize',2)
plot([0 Qi1(1)], [0 Qi1(2)], 'r','linewidth',2)
plot([0 Qi2(1)], [0 Qi2(2)], 'g','linewidth',2)
plot([0 Qi1(1)], [0 Qi1(2)], 'k--','linewidth',2)
plot([0 Qi2(1)], [0 Qi2(2)], 'k--','linewidth',2)
plot(-[0 Qi1(1)], -[0 Qi1(2)], 'r','linewidth',2)
plot(-[0 Qi2(1)], -[0 Qi2(2)], 'g','linewidth',2)
plot(-[0 Qi1(1)], -[0 Qi1(2)], 'k--','linewidth',2)
plot(-[0 Qi2(1)], -[0 Qi2(2)], 'k--','linewidth',2)
title(['separation index = ' num2str(sum([QSubspaces.costFn]))])
end
%%