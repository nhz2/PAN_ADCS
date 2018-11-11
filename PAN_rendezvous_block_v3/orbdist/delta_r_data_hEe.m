function [deltars,deltarns] = delta_r_data_hEe(x_t,x_f,muE,nsteps)
%delta_r_data_hEe - compute delta_r vectors and norm(delta_r) values using
%hEe method

h_t = x_t(1:3);
E_t = x_t(4);
ev_t = x_t(5:7);
a_t = -muE/(2*E_t);
e_t = norm(ev_t);
hh_t = h_t/norm(h_t);
eh_t = ev_t/norm(ev_t);
orb_t = [a_t;e_t;hh_t;eh_t];

h_f = x_f(1:3);
E_f = x_f(4);
ev_f = x_f(5:7);
a_f = -muE/(2*E_f);
e_f = norm(ev_f);
hh_f = h_f/norm(h_f);
eh_f = ev_f/norm(ev_f);
orb_f = [a_f;e_f;hh_f;eh_f];

d_orb = orb_t - orb_f;
deltars = zeros(3,nsteps);
deltarns = zeros(1,nsteps);
nuspan = linspace(0,2*pi,nsteps);

for i=1:nsteps
    deltars(:,i) = delta_r_hEe(orb_t,d_orb,nuspan(i));
    deltarns(i) = norm(deltars(:,i));
end
end