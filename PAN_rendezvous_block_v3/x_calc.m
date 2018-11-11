function [x] = x_calc(r,v,muE)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
h=cross(r,v);
E=norm(v)^2/2-muE/norm(r);
e=(cross(v,h)/muE-r/norm(r));
x=[h(1);h(2);h(3);E;e(1);e(2);e(3)];
end

