function dr_maxs = delta_r_maxs(x_t,x_f,muE,nsteps)
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
dr_maxs = zeros(4,nsteps);
nuspan = linspace(0,2*pi,nsteps);

for i=1:nsteps
    deltar_i = delta_r_hEe(orb_t,d_orb,nuspan(i));
    
    for j=1:3
        if (abs(deltar_i(j)) > abs(dr_maxs(j)))
            dr_maxs(j) = deltar_i(j);
        end
    end
    
    if (norm(deltar_i) > dr_maxs(4))
        dr_maxs(4) = norm(deltar_i);
    end
end
end