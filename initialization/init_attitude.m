function init= init_attitude(init)
% RWA Properties
init.FWR.h= [1e-3 1e-3 1e-3]';    % initial angular velocity of the wheels
%init.LDR.h= [1.125e-3;1.125e-3;1.125e-3]*1.5;    % initial angular velocity of the wheels
%init.h_des= [1.125e-3;1.125e-3;1.125e-3];
init.LDR.h= [1e-3 1e-3 1e-3]';

% For Time-Domain Solution
%init.wbn=pi/30*[30 15 1]';
%init.FWR.wbn=pi/30*[0 0 0]';
init.FWR.wbn=[0 0 -0.001]';
%init.FWR.w_ref=pi/30*[0 0 0];
init.FWR.q_init=[0 0 0 1]';
%init.FWR.q_ref=[0 0 0 1]';

init.LDR.wbn=pi/100*[1 0 0]';
init.LDR.w_ref=pi/30*[0 0 0];
%init.LDR.q_init=[0 0 0 1]';
init.LDR.q_ref=[0 0 0 1]';
a= [1; 0; 0];
th= pi/180 * 91;
%init.LDR.q_init=[cos(pi/175/2)*sin(pi/175/2), -sin(pi/175/2)*sin(pi/175/2), 0,sin(pi/175/2)];
init.LDR.q_init=[a*sin(th/2);cos(th/2)];
init.FWR.q_init=[a*sin(th/2);cos(th/2)];
v1= abs(init.FWR.P);
v2= abs(cross(init.FWR.V,init.FWR.P));
v3= abs(cross(v1,v2));
init.DCM_be= [v1,v2,v3];
%

%{
% For Time-Domain Solution
wbn_init=pi/30*[0 0 0]';
w_ref=pi/30*[0 0 0];
q_init=[0 0 0 1]';
q_ref=[0 0 0 1]';
%}