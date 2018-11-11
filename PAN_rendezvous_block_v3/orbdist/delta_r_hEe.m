function delta_r = delta_r_hEe(orb_ref,dorb,nu)
%delta_r_hEe - find change in position at a given true anomaly from reference
%orbit for change in state dx

a = orb_ref(1);
e = orb_ref(2);
hh = orb_ref(3:5);
eh = orb_ref(6:8);
qh = cross(hh,eh);

da = dorb(1);
de = dorb(2);
dhh = dorb(3:5);
deh = dorb(6:8);

Q = [eh qh hh];
Prh = [cos(nu); sin(nu); 0];
Irh = Q*Prh;

dr_a = (1-e^2)/(1+e*cos(nu))*Irh*da;
dr_e = (-2*a*e*(1+e*cos(nu))-a*(1-e^2)*cos(nu))/(1+e*cos(nu))*Irh*de;
dr_hh = a*(1-e^2)/(1+e*cos(nu))*(sin(nu)*(crs(eh)+eye(3)))*dhh;
dr_eh = a*(1-e^2)/(1+e*cos(nu))*(cos(nu)*eye(3)-sin(nu)*crs(hh))*deh;

delta_r = dr_a+dr_e+dr_hh+dr_eh;

end

