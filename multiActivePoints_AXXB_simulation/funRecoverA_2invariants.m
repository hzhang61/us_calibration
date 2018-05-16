%% Definition
% Input: 
%   v: initial values
%   subA: subscript As containing 1DoF uncertainty
%   subB: subscript Bs
%   subA1: fully recontructed subscript A for the first pose
%
% Output: minimization cost function to find full subscript As using two
% invariants

% /***************************************************************************
% Copyright 
% MUSiiC Laboratory
% Haichong Zhang,Emad M Boctor
% Johns Hopkins University
% 
% For commercial use/licensing, please contact Emad Boctor, Ph.D. at eboctor@jhmi.edu.
% ***************************************************************************/

function F = funRecoverA_2invariences(v, subA, subB, subA1)

% initialization
m = size(subA,3);
subA(:,:,1) = subA1;
v(1) = 0;
RotA = zeros(4,4,m);
RotA(:,:,1) = buildT([0 0 0 0 0 0]);
nsubA(:,:,1) = (inv(subA(:,:,1))*RotA(:,:,1));

for i = 1:m
    RotA(:,:,i) = buildT([0 0 v(i) 0 0 0]);
    nsubA(:,:,i) = (inv(subA(:,:,i))*RotA(:,:,i));
end

% computing superscript As and Bs
count = 0;
for i = 1:m
    for j = 1:m
        if i ~= j
            count = count + 1;
            supA(:,:,count) = inv(nsubA(:,:,i))*(nsubA(:,:,j));
            supB(:,:,count) = inv(subB(:,:,i))*subB(:,:,j);
        end
    end
end

% minimizing invariants
count = 0;
for i = 1:m
    for j = 1:m
        if i ~= j
            count = count + 1;
            [theta_B(count), N_B(:,:,count), d_B(count),p_B(:,count)]=param_extract2(supB(:,:,count));
            n_B(1,count) = N_B(3,2,count);
            n_B(2,count) = N_B(1,3,count);
            n_B(3,count) = N_B(2,1,count);
            temp = supA(:,:,count);
            [theta_A(count), N_A(:,:,count), d_A(count),p_A(:,count)]=param_extract2(temp);
            n_A(1,count) = N_A(3,2,count);
            n_A(2,count) = N_A(1,3,count);
            n_A(3,count) = N_A(2,1,count);
            F(count*2-1) = (theta_B(count)-theta_A(count));
            F(count*2) = (d_B(count)-d_A(count));
        end
    end
end