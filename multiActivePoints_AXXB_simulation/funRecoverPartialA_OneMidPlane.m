%% Definition
% Input:
%   v: initial values
%   model: phantom model
%   sp: segmented point position
%   errorOOP: out of plane error to be simulated
%
% Output: minimization cost function to find full subscript A

% /***************************************************************************
% Copyright 
% MUSiiC Laboratory
% Haichong Zhang,Emad M Boctor
% Johns Hopkins University
% 
% For commercial use/licensing, please contact Emad Boctor, Ph.D. at eboctor@jhmi.edu.
% ***************************************************************************/

%% Function
function F = funRecoverPartialA_OneMidPlane(v, model, sp, errorOOP)

A = buildT(v(1:6));

n = size(model,2);
nsp = zeros(3,n);
count = 0;
for point = 1:n
    count = count + 1;
    if point == 1
        v(6+point) = errorOOP;
    end
    nsp(1,point) = sp(1,point);
    nsp(2,point) = sqrt(sp(2,point)^2-v(6+point)^2);
    nsp(3,point) = v(6+point);
    cp(:,point) = (A)*model(:,point);
    
    F(count) = norm(nsp(1,point)-cp(1,point));
    count = count + 1;
    F(count) = norm(nsp(2,point)-cp(2,point));
    count = count + 1;
    F(count) = norm((nsp(3,point))-(cp(3,point)));
end