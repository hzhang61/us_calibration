%% Definition
% Input: 3 by 3 rotational transformation 
%
% Output: Euler angle representation of a rotation matrix

% /***************************************************************************
% Copyright 
% MUSiiC Laboratory
% Haichong Zhang,Emad M Boctor
% Johns Hopkins University
% 
% For commercial use/licensing, please contact Emad Boctor, Ph.D. at eboctor@jhmi.edu.
% ***************************************************************************/

function [x,y,z] = decompose_rotation_d(R)
	x = (atan2(R(3,2), R(3,3))/pi)*180;
	y = (atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)))/pi)*180;
	z = (atan2(R(2,1), R(1,1))/pi)*180;
end