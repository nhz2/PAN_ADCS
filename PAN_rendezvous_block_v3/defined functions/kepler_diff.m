% Matthew Walsh
% Kaplerian Orbit Differential Equations
% Differential Equation Function

function [ zdot ] = kepler_diff( t,z,mu,u )
%Kepler_diff - computes rate of change for state vector z according to 
% differential equations of gravitation about a central body with 
% additional forcesfor a small satellite, central body is assumed to be 
% stationary, mass of satellite is assumed to be negligible
% z contains the positions and velocities for satellite
% mu is the gravitational parameter
% u is the vector of other forces acting on the satellite divided by its 
%    mass [ux uy uz]'
% z should contain states 
%    [x y z xd yd zd]'

a=z(4);
b=z(5);
c=z(6);

x=z(1);
y=z(2);
z=z(3);

r=sqrt(x^2+y^2+z^2);

zdot=zeros(6,1);
zdot(1)=a;      % xdot
zdot(2)=b;      % ydot
zdot(3)=c;      % zdot
zdot(4)=-mu*x/r^3+u(1); % x1doubledot
zdot(5)=-mu*y/r^3+u(2); % x2doubledot
zdot(6)=-mu*z/r^3+u(3); % y1doubledot


end