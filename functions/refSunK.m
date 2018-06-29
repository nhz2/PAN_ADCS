function S= refSunK(t,tp,yearT,e,q)%#codegen
%refSun Return the normalized vector from earth to sun in ECI
%       A very simple model of the earth orbiting around the sun in a keplerian orbit.
%   t(positive double): delta time in seconds
%   tp(double): delta time when the earth is at perihelion in
%       seconds
%   yearT(positive double): earths orbital period in seconds
%   e(positive double): eccentricity of earths orbit
%   q(quaterinion): quaternion to rotate from Sun to ICRF
%TODO use custom quat library
%   
%   
%find M in rads
M= 2*pi*(t-tp)/yearT;
% find E in rads using fixed point iteration see
% https://en.wikipedia.org/wiki/Kepler%27s_equation
E= M;
E= M+e*sin(E);
% for n= 0:100
%     E= M+e*sin(E);
% end
%get xses and yses
xses= cos(E)-e;
yses= sqrt(1-e*e)*sin(E);
d= sqrt(xses*xses+yses*yses);
xses= xses/d;
yses= yses/d;

%rotate to ICRS aka ECI
qu= quaternion(q(4),-q(1),-q(2),-q(3));
%S= zeros([3,1]);
S= -rotatepoint(qu,[xses,yses,0]).';
end

