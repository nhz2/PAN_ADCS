function S= refSun(t,angle0,tilt,yearT)%#codegen
%refSun Return the normalized vector from earth to sun in ECI
%       A very simple model of the earth rotating around the sun with a
%       fixed axial tilt.
%   t(positive double): delta time in seconds
%   
%   
S=[cos(2*pi*t/yearT+angle0)+0.015*sin(4*pi*t/yearT)+0.01*cos(2*pi*t/yearT)-0.01*cos(4*pi*t/yearT)-0.005; sin(2*pi*t/yearT+angle0)*cos(tilt); sin(2*pi*t/yearT+angle0)*sin(tilt)];
    


end

