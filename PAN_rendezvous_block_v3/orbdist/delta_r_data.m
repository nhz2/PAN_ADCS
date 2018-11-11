function [deltars,deltarns] = delta_r_data(x_t,x_f,muE,nsteps)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[a_t, e_t, I_t, Om_t, w_t] = x2orb(x_t,muE);
els_t = [a_t; e_t; I_t; Om_t; w_t];
[a_f, e_f, I_f, Om_f, w_f] = x2orb(x_f,muE);
els_f = [a_f; e_f; I_f; Om_f; w_f];
d_els = els_t - els_f;
deltars = zeros(3,nsteps);
deltarns = zeros(1,nsteps);
nuspan = linspace(0,2*pi,nsteps);
 
for i=1:nsteps
    deltars(:,i) = delta_r(els_t,d_els,nuspan(i));
    deltarns(i) = norm(deltars(:,i));
end
end

