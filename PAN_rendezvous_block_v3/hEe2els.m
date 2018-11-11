function [a,e,i,Om,w] = hEe2els(hEe)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    h=hEe(1:3);
    E=hEe(4);
    ev=hEe(5:7);
    
    a=-muE/(2*E);
    e=norm(ev);
    i=acos(h(3)/norm(h0_lf));
    
    n=cross([0;0;1],h);
    Om=acos(n(1)/norm(n));
    q34=n(2)<0;
    Om(q34)=2*pi-Om(q34);
    
    if (e0v_f(3)<0)
        w=2*pi-acos(dot(n,ev)/(norm(n)*norm(ev)));
    else
        w=acos(dot(n,ev)/(norm(n)*norm(ev)));
    end
end

