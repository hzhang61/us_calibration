%% Definition
% Input: the scale of the phantom model
%
% Output: phantom model with 12 points

% /***************************************************************************
% Copyright 
% MUSiiC Laboratory
% Haichong Zhang,Emad M Boctor
% Johns Hopkins University
% 
% For commercial use/licensing, please contact Emad Boctor, Ph.D. at eboctor@jhmi.edu.
% ***************************************************************************/

function [model] = funPhantomModel(m)

scale = m;
tmodel(1,:) = [0 0 0 1];
tmodel(2,:) = [5 5 0 1];
tmodel(3,:) = [0 5 0 1];
tmodel(4,:) = [5 0 0 1];
tmodel(5,:) = [0 0 5 1];
tmodel(6,:) = [5 5 5 1];
tmodel(7,:) = [0 5 5 1];
tmodel(8,:) = [5 0 5 1];
tmodel(9,:) = [0 0 10 1];
tmodel(10,:) = [5 5 10 1];
tmodel(11,:) = [0 5 10 1];
tmodel(12,:) = [5 0 10 1];
tmodel(:,1:3) = tmodel(:,1:3)*scale;
model = tmodel';

end
