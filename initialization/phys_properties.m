function init= phys_properties(init)

% Earth
init.Earth.mu= 3.986004418e14;   % Gravitational Constant [m^3/s^2]
init.Earth.M = 5.97237e24;       % Mass [kg]
init.Earth.R = 6371e3;           % Mean radius [m]

% Sun
init.Sun.M = 1.98855e30;         % Mass [kg]
init.Sun.R = 695700e3;           % Radius [m]

% 3U+ CubeSat
init.CubeSat.Dim = [10 30 10]*0.01;  % [m]
init.CubeSat.M = 4; % [kg]
m= init.CubeSat.M;
init.I=[1/12*m*(0.3^2+0.1^2) 0 0
        0 1/12*m*(0.3^2+0.1^2) 0
        0 0 1/12*m*(0.1^2+0.1^2)];
init.invI=inv(init.I);
init.c=100;

%Aerodynamic/Solar Properties
init.Cd= 2.0;
%drag coefficient
%faces +x,-x,+y,-y,+z,-z 
init.areas= [0.03 0.03 0.03 0.03 0.01 0.01];
%areas([1, N]  vector): areas of the planes, units m^2
init.unitnormals= [ 1 -1 0  0 0  0
                    0  0 1 -1 0  0
                    0  0 0  0 1 -1];
%unitnormals([3, N]  vectors): normalzed unit vectors of planes
%   these point out of the plane into space, in body frame
init.comtocops= [ 5 -5 0  0  0   0
                  0  0 5 -5  0   0
                  0  0 0  0 10 -10]*0.01;
%comtocops([3, N]  vectors): vector distance from center of
%   mass to center of pressure of each plane, units m, in body frame
init.rss= [.4 .4 .4 .4 .4 .4];
%rss([1, N]  vector): surface specular reflectance coefficients.
init.rds= [.4 .4 .4 .4 .4 .4];
%rds([1, N]  vector): surface diffuse reflectance coefficients.



% Wheel Properties (Maxon EC 45 Flat 50W 24V)
init.ratemax= 677.0;      % Maximum speed rad/s
init.rampmax= 304.5;        %Maximum ramp in rad/s/s
init.wheel_I= 135*1e-7;        % Rotor Inertia kg*m^2
init.ratedes= 50.0;     %Desired wheel speed (rad/s)

%Magnetic torquer
init.m_mt= 0.08; %Am^2 Maximum moment of a single rod

