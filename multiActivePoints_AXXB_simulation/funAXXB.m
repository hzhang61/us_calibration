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
function F = funAXXB(v, A, B)

X = buildT(v);
m = size(A,3);
for j = 1:m
    temp1 = (A(:,:,j))*X;
    temp2 = X*B(:,:,j);
    F(6*j-5:6*j-3) = decompose_rotation_d(temp1(1:3,1:3)*inv(temp2(1:3,1:3)));
    F(6*j-2:6*j) = temp1(1:3,4)-temp2(1:3,4);
end