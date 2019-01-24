% PAN rendezvous GNC block test Master Script
% Version 2
% Matt Walsh

% global constants
RE = 6378137; % radius of earth in m
muE = 3.98600441500E+14; % gravitational parameter of earth m^3/s^2
J2 = 1.08263566655E-03; % oblateness constant

%% Initialization

% number of steps for optimal dV calculation
nsteps = 10000;
nsteps_dr = 1000;

% maneuver switch; set to 1
manswitch = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Specification                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System Properties
% spacecraft mass in kg
m_follower = 5;

% hardware limits
% subsystems not fully implemented so parameters given are placeholders

% thrust range in N
Thrust_min=4.9e-3;
Thrust_max=5.8e-3;

% ontime range in s
ontime_min = 0.010;
ontime_max = 1;

% impulse per pulse range for thrusters in Ns
% Imp_range=[Thrust_min*ontime_min Thrust_max*ontime_max];
Imp_range = [4.9e-5 5.8e-3];     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time parameters                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_update=1200; % seconds between th_opt updates
t_wait=7200; % initial wait time - time after activation maneuvers begin
man_time=1200; % length of maneuver calculation and execution cycle (s)
t_end=3600*24*20; % time to end simulation
T_man = 1200; % time between maneuvers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tolerances                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Distance for close-proximity rendezvous
d_target=2000; % m; distance for docking range

% tolerances for primary orbit matching
rtols = [500;500;50];
rntol = 500;

%% Orbit Specification

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbital Elements                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% els_t = [RE+405000; 4000/(2*(RE+405000)); 51.6403*pi/180; 189.1906*pi/180;
%         283.7256*pi/180; pi/3];
% els_f = [RE+403000; .0005; 51.65*pi/180; 185*pi/180; 283.72*pi/180; 2*pi/3];
% 
% a0_t=els_t(1); % semimajor axis
% e0_t=els_t(2); % eccentricity
% i0_t=els_t(3); % inclination
% w0_t=els_t(4); % argument of periapsis
% Om0_t=els_t(5); % longitude of ascending node
% th0_t=0; % true anomaly
% 
% a0_f=els_f(1); % semimajor axis
% e0_f=els_f(2); % eccentricity
% i0_f=els_f(3); % inclination
% w0_f=els_f(4); % argument of periapsis
% Om0_f=els_f(5); % longitude of ascending node
% th0_f=0; % true anomaly

% follower
r0_f=[6252006.47153675;-2568969.80746986;407208.247032049];
v0_f=[1479.79410947239;4550.20536946861;5999.73074766994];

% target
r0_t=[6269648.63827392;-814926.681168529;2424333.57589868];
v0_t=[-1379.44461861257;5322.84133137884;5351.02269477794];

%% Simulation setup and drift phase

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% orbit setup for sim                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% follower
% [r0_f,v0_f]=orb2rv_s(a0_f*(1-e0_f^2),e0_f,i0_f,Om0_f,w0_f,th0_f,muE);
% h0_f=cross(r0_f,v0_f);
% E0_f=norm(v0_f)^2/2-muE/norm(r0_f);
% e0v_f=cross(v0_f,h0_f)/muE-r0_f/norm(r0_f);
% x0_f=[h0_f;E0_f;e0v_f];

h0_f=cross(r0_f,v0_f);
E0_f=norm(v0_f)^2/2-muE/norm(r0_f);
ev0_f=cross(v0_f,h0_f)/muE-r0_f/norm(r0_f);
x0_f=[h0_f;E0_f;ev0_f];
[a0_f,e0_f,i0_f,Om0_f,w0_f,th0_f]=rv2orb(r0_f,v0_f,muE);

% target
% [r0_t,v0_t]=orb2rv_s(a0_t*(1-e0_t^2),e0_t,i0_t,Om0_t,w0_t,th0_t,muE);
% h0_t=cross(r0_t,v0_t);
% E0_t=norm(v0_t)^2/2-muE/norm(r0_t);
% e0v_t=cross(v0_t,h0_t)/muE-r0_t/norm(r0_t);
% x0_t=[h0_t;E0_t;e0v_t];

h0_t=cross(r0_t,v0_t);
E0_t=norm(v0_t)^2/2-muE/norm(r0_t);
ev0_t=cross(v0_t,h0_t)/muE-r0_t/norm(r0_t);
x0_t=[h0_t;E0_t;ev0_t];
[a0_t,e0_t,i0_t,Om0_t,w0_t,th0_t]=rv2orb(r0_t,v0_t,muE);

%% Pre simulation calculations

% Weightings for optimal thrust calculation

K_F=diag([1/norm(h0_t);1/norm(h0_t);1/norm(h0_t);1/abs(E0_t);1;1;1]);
L_F=eye(3);
% K_F=diag([0;0;0;1;0;0;0]);

% Assign impulse limits
Imp_min = Imp_range(1);
Imp_max = Imp_range(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SimIn.Id.Name = 'hEe_orbit_rendezvous';
SimIn.Id.Version = '3';
SimIn.OrbMatch.EOM_IC.r0_t = r0_t;
SimIn.OrbMatch.EOM_IC.r0_f = r0_f;
SimIn.OrbMatch.EOM_IC.v0_t = v0_t;
SimIn.OrbMatch.EOM_IC.v0_f = v0_f;
SimIn.Resources.dV_used = 0;

