%% Definition
%function to extract Plucker parameters from a g \in SE(3)

%% Function
function [theta, N, d, p]=param_extract2(X)

%extract theta---------------------------------
theta=acos((trace(X(1:3,1:3))-1)/2);
%----------------------------------------------

%extract N---------------------------------
if theta<7e-2
    N=so3_vec(X(1:3,4)/norm(X(1:3,4)));
else
    N=(X(1:3,1:3)-X(1:3,1:3)')/(2*sin(theta));
end
%----------------------------------------------

%extract d---------------------------------
d=dot(X(1:3,4),so3_vec(N));
%----------------------------------------------

%extract p---------------------------------
n=so3_vec(N);
u=1/sqrt(n(1)^2+n(2)^2)*[-n(2); n(1); 0];
c=inv([1-cos(theta) sin(theta); -sin(theta) 1-cos(theta)])*([dot(X(1:3,4),u); dot(X(1:3,4),(cross(n,u)))]);
p=c(1)*u+c(2)*cross(n,u);
%----------------------------------------------
end