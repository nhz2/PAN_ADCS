function S= refSunK(t,tp,yearT,e,tilt,anomeq,N)%#codegen
%refSun Return the normalized vector from earth to sun in ECI
%       A very simple model of the earth orbiting around the sun in a keplerian orbit.
%   t(positive double): delta time in seconds
%   tp(double): delta time when the earth is at perihelion in
%       seconds
%   yearT(positive double): earths orbital period in seconds
%   e(positive double): eccentricity of earths orbit
%   tilt(double): delta time in seconds
%   anomeq(double): true anomoly of earth at equoniox, when sun to earth
%       vector is parrallel to  (1, 0, 0)ICRS units rad
%   N(positive int): number of iterations in inverse kepler equation
%   
%   
%   
%find M in rads
M= 2*pi*(t-tp)/yearT;
% find E in rads using fixed point iteration see
% https://en.wikipedia.org/wiki/Kepler%27s_equation
E= M;
for n= 0:N
    E= M+e*sin(E);
end
%get xses and yses
xses= cos(E)-e;
yses= sqrt(1-e*e)*sin(E);
d= sqrt(xses*xses+yses*yses);
xses= xses/d;
yses= yses/d;

%rotate to ICRS aka ECI
cp= cos(anomeq);
sp= sin(anomeq);
ct= cos(tilt);
st= sin(tilt);
A= [ cp      sp     0;
    -ct*sp   ct*cp  st;
     st*sp  -st*cp  ct];
S= (A*[xses; yses; 0]);
end

