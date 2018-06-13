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
init.JD= (1461 * (Y + 4800 + (M - 14)/12))/4 +(367 * (M - 2 - 12 * ((M - 14)/12)))/12 - (3 * ((Y + 4900 + (M - 14)/12)/100))/4 + D - 32075;
init.GPSTime= 1213056018;
%seconds, using https://losc.ligo.org/gps/

dataIn= double('1201,20304!1202,302350!1301,12352!');
dataArray= [12 32 2124.6 1234 12314 5 6 5 4 3 6 7 8 9 23];