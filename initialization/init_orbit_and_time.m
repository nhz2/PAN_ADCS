function init= init_orbit_and_time(init)
a= 6760636.6;   % semimajor axis
e= 0;           % eccentricity
p= a*(1-e);     % semilatus rectum
i_LDR= 0.01*pi/180;   % inclination angle
i_FWR= 0*pi/180;   % inclination angle
O=0;            % right ascension of the ascending node (longitude)
o=0;            % Argument of Perigee
nu_LDR=5*pi/180;           % True anamoly
nu_FWR=0*pi/180;           % True anamoly
[r_LDR, v_LDR]= orb2rv_s(p,e,i_LDR,O,o,nu_LDR,init.Earth.mu);
[r_FWR, v_FWR]= orb2rv_s(p,e,i_FWR,O,o,nu_FWR,init.Earth.mu);

% CubeSat Orbital Initial conditions
init.time= 0;
%init.LDR.P = [6760636.6 0 0]';   % Px Py Pz [m]
%init.LDR.V = [0 7.678477071156148e+03 0]';     % Vx Vy Vz
init.LDR.P = r_LDR;
init.LDR.V = v_LDR;

init.FWR.P = r_FWR;   % Px Py Pz [m]
init.FWR.V = v_FWR;     % Vx Vy Vz

init.errorflag = 'Error';
init.deltaT = 'sec';

% Start Date UTC
init.year= 2018;
init.month= 6;
init.day= 15;
init.hour= 0;
init.min= 0;
init.sec= 0;
Y= init.year;
M= init.month;
D= init.day;
init.JD= juliandate(Y,M,D,init.hour,init.min,init.sec);
init.GPSTime= 1213056018;
init.decyear= decyear(Y,M,D);
%seconds, using https://losc.ligo.org/gps/

%Sun earth orbit init constants 
%https://en.wikipedia.org/wiki/Earth%27s_orbit
init.S_axialtilt= 23.43689*pi/180.0; 
%radians see https://en.wikipedia.org/wiki/Axial_tilt
init.S_eqangle= 1.459420149;
%init angle of earth from equaniox 
init.S_T= 365.256363004*86400.0;
%earths orbital period units s
init.S_e= 0.0167086;
%earths orbital eccentricity
init.S_tp= (juliandate(2018,1,3,5,35,0)- init.JD)*86400.0;
%delta time when earth is at perihelion units s
%http://aa.usno.navy.mil/seasons?year=2018&tz=+0
r= planetEphemeris(juliandate(2018,1,3,5,35,0),'Sun','Earth','421','AU');
r= r/sqrt(dot(r,r));
init.S_anomeq= acos(r(1));
%true anomoly of earth at equoniox, when sun to earth
%vector is parrallel to  (1, 0, 0)ICRS units rad


dataIn= double('1201,20304!1202,302350!1301,12352!');
dataArray= [12 32 2124.6 1234 12314 5 6 5 4 3 6 7 8 9 23];