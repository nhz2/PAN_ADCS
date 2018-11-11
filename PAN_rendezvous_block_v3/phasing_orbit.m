function [r_tgt,v_tgt,x_tgt] = phasing_orbit(els_current,a_step,muE)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

a_current=els_current(1);
e_current=els_current(2);
i_current=els_current(3);
Om_current=els_current(4);
w_current=els_current(5);
th_current=els_current(6);

r1_tgt = a_current*(1-e_current);
r2_tgt = 2*a_step - r1_tgt;
e_tgt = abs((r1_tgt-r2_tgt)/(2*a_step));

w_tgt = w_current;
th_tgt = th_current;

if (r1_tgt > r2_tgt)
    w_tgt = w_tgt + pi;
    th_tgt = th_tgt - pi;
    if (w_tgt > 2*pi)
        w_tgt = w_tgt - 2*pi;
    end
    
    if (th_tgt < 0)
        th_tgt = th_tgt + 2*pi;
    end
end

[r_tgt,v_tgt]=orb2rv_s(a_step*(1-e_tgt^2),e_tgt,i_current,Om_current,w_tgt,th_tgt,muE);
h_tgt = cross(r_tgt,v_tgt);
E_tgt = -muE/(2*a_step);
ev_tgt = cross(v_tgt,h_tgt)/muE - r_tgt/norm(r_tgt);
x_tgt = [h_tgt;E_tgt;ev_tgt];
end

