init= struct();
init= phys_properties(init);
init= init_orbit_and_time(init);
init= init_attitude(init);
init= sensor_prop(init);
init= driver_prop(init);

%GNC Initializations
muE= init.Earth.mu; %gravitational parameter for the earth
m_follower= init.CubeSat.M; %mass of the spacecraft
Thrust_max= 1; %maximum thrust per impulse
Thrust_min= 1; %minimum thrust per impulse
ontime_max= 1; %maximum ontime
ontime_min= 1; %minimum ontime
Imp_max= 1; %maximum impulse per maneuver
Imp_min= 1; %minimum impulse per maneuver
%weighting matrices for optimization calculation
h0_t = 1;%angular momentum
E0_t = 1;%Energy
K_F= diag([1/norm(h0_t);1/norm(h0_t);1/norm(h0_t);1/abs(E0_t);1;1;1]);

L_F= eye(3);
match_tol=0.01; %tolerance for orbit matching, how the block determines which phase of the rendezvous it’s in


% These gains are designed for infinitesimal errors.  The high rate at t=0
% causes a nonlinear response that nonetheless settles (per Lyapunov)
Kp=7.5e-7;     
Kd=75e-7;    
Kp=7.5e-4;     
Kd=75e-4; 
