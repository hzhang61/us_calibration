%% Definition
% Input:
%   A: superscript As
%   B: superscript Bs
%
% Output: minimization cost function to find X

% /***************************************************************************
% Copyright 
% MUSiiC Laboratory
% Haichong Zhang,Emad M Boctor
% Johns Hopkins University
% 
% For commercial use/licensing, please contact Emad Boctor, Ph.D. at eboctor@jhmi.edu.
% ***************************************************************************/

%% Function
function [ X, residual ] = solveAXXB( A, B)
minCost = Inf;
opts=  optimset('display','off','MaxFunEvals',5000,'MaxIter',1000);

for i = 1:5
    initV = [2*pi*rand(6,1)];
    [v,residual] = lsqnonlin(@funAXXB,initV,[],[],opts,A,B);
    if residual < minCost
        minCost = residual;
        finalV = v;
    end
end
X = buildT(finalV(1:6));
end