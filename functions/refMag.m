function B= refMag(t,x)%#codegen
%refMag Return the magnetic feild vector in ECEF units gauss
%   t(positive double): decyear
%   x([3, 1]  vector): position in ECEF coords, units meters
%   
%   using model described in WMM2015_Report
%       Chulliat, A., S. Macmillan, P. Alken, C. Beggan, M. Nair, 
%       B. Hamilton, A. Woods, V. Ridley, S. Maus and A. Thomson, 2015, 
%       The US/UK World Magnetic Model for 2015-2020: Technical Report, 
%       National Geophysical Data Center, NOAA. doi: 10.7289/V5TB14V7

B=zeros([3,1]);


coder.updateBuildInfo('addSourceFiles','GeoMag.c');
coder.cinclude('GeoMag.h');

xstruct= struct('x',x(1),'y',x(2),'z',x(3));
coder.cstructname(xstruct,'MAGtype_CoordECEF','extern','HeaderFile','GeoMag.h');
rstruct= struct('Bx',B(1),'By',B(2),'Bz',B(3));
coder.cstructname(rstruct,'MAGtype_MagneticResults','extern','HeaderFile','GeoMag.h');
coder.ceval('MAG_geomag',xstruct,t,coder.ref(rstruct));
B(1)= rstruct.Bx;
B(2)= rstruct.By;
B(3)= rstruct.Bz;

end





