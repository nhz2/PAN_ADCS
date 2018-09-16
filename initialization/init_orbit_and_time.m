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
init.LDR.H = cross(r_LDR,v_LDR); %angular momuntum/mass in ECI units m^2/s
%https://en.wikipedia.org/wiki/Angular_momentum


init.FWR.P = r_FWR;   % Px Py Pz [m]
init.FWR.V = v_FWR;     % Vx Vy Vz

init.errorflag = 'Error';
init.deltaT = 'sec';

% Start Date UTC, don't change this without changing the ecef to eci blocks
init.year= 2017;
init.month= 10;
init.day= 20;
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
init.S_T= 365.256363004*86400.0;
%earths orbital period units s
init.S_e= 0.0167086;
%earths orbital eccentricity
init.S_tp= (juliandate(2018,1,3,5,35,0)- init.JD)*86400.0;
%delta time when earth is at perihelion units s
%http://aa.usno.navy.mil/seasons?year=2018&tz=+0
x= planetEphemeris(juliandate(2018,1,3,5,35,0),'Sun','Earth','421','AU');
x= x/sqrt(dot(x,x));
p= planetEphemeris(juliandate(2018,3,3,5,35,0),'Sun','Earth','421','AU');
p= p/sqrt(dot(p,p));
z= cross(x,p);
z= z/sqrt(dot(z,z));
init.S_q= TRIAD([0; 0; 1],[1; 0; 0],z.',x.');
%quaternion to rotate from Sun to ICRF


dataIn= double('1201,20304!1202,302350!1301,12352!');
dataArray= [12 32 2124.6 1234 12314 5 6 5 4 3 6 7 8 9 23];