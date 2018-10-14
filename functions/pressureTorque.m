function [force,torque] = pressureTorque(pressure,unitvelocity,areas,unitnormals,comtocops,rss,rds)
%PRESSURETORQUE Return the force and torque of a pressure on some faces
%   See page 324, table 10-18 in SMAD 
%       and page 108, equation 3.163 in Fundamentals of Spacecraft ADC
%       for equations this is based on.
%   The force on a face is applied at its center of pressure(cop).
%   The force is zero if unitvelocity dot unitnormal is nonnegative.
%   Other wise the force is (unitvelocity dot unitnormal)*pressure*area*
%      ((2*rs+2/3*rd)*unitnormal+(rs-1)*unitvelocity)
%   Args:
%       pressure(double): the pressure constant, units Pa
%           For solar pressure, 4.644E-6 Pa
%           For air, 0.5*rho*Cd*v^2
%       unitvelocity([3, 1]  vector): normalzed velocity of particles
%           For solar pressure, vector from sun to sat
%           For air, negative sat velocity relative to the air
%       areas([1, N]  vector): areas of the planes, units m^2
%       unitnormals([3, N]  vectors): normalzed unit vectors of planes
%           these point out of the plane into space.
%       comtocops([3, N]  vectors): vector distance from center of
%           mass to center of pressure of each plane, units m
%       rss([1, N]  vector): surface specular reflectance coefficients.
%           0 for air.
%       rds([1, N]  vector): surface diffuse reflectance coefficients.
%           0 for air.
N= length(areas);
assert(isequal(size(unitvelocity),[3 1]),"bad dimesions");
assert(isequal(size(comtocops),size(unitnormals),[3 N]),"bad dimesions");
assert(isequal(size(areas),size(rss),size(rds),[1 N]),"bad dimesions");
force= [ 0; 0; 0;];
torque= [ 0; 0; 0;];
for i = 1:N
    d= dot(unitnormals(i),unitvelocity);
    if d>=0
        continue;
    end
    f= d*pressure*areas(i)*((2*rss(i)+2/3*rds(i))*unitnormals(i)+(rss(i)-1)*unitvelocity);
    force= force+f;
    torque= torque+ cross(comtocops(i),f);  
end
end