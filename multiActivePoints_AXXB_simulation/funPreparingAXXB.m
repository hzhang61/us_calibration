%% Definition
% Input: number of poses to be generated. X and Y are arbiturary
% determined
%
% Output: Randomly generated A and B

% Variables
% subA: subscript A
% supA: superscript A
% subB: subscript B
% supB: superscript B
% X: transformation from tracker frame to US image frame
% Y: transformation from tracker base frame to phantom frame

% /***************************************************************************
% Copyright 
% MUSiiC Laboratory
% Haichong Zhang,Emad M Boctor
% Johns Hopkins University
% 
% For commercial use/licensing, please contact Emad Boctor, Ph.D. at eboctor@jhmi.edu.
% ***************************************************************************/

%% Function
function [subA, supA, subB, supB, X, Y] = funPreparingAXXB(poses)

% generating randomized subscript As
for i = 1:poses
    subA(:,:,i) = buildT([40*randn(1,3) abs(10*randn(1,1)) abs(20*randn(1,1))+100 0]);
end

% arbitrary determined X and Y
X = buildT([29 41 14 14 53 109]);
Y = buildT([-52 56 35 245 434 254]);

% subscript Bs are determined based on A, X, and Y
for i = 1:poses
    subB(:,:,i) = Y*inv(subA(:,:,i))*(X);
end

% superscript As and Bs are generated from subscript As and Bs
count = 0;
for i = 1:poses
    for j = 1:poses
        if i ~= j
            count = count + 1;
            supA(:,:,count) = (subA(:,:,i))*inv(subA(:,:,j));
            supB(:,:,count) = inv(subB(:,:,i))*subB(:,:,j);
        end
    end
end