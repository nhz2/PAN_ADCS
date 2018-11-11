function delta_r = delta_r(els_ref,d_els,nu)
%delta_r - find change in position at a given true anomaly from reference
%orbit for change in orbital elements d_els

a = els_ref(1);
e = els_ref(2);
I = els_ref(3);
O = els_ref(4);
w = els_ref(5);

da = d_els(1);
de = d_els(2);
dI = d_els(3);
dO = d_els(4);
dw = d_els(5);

r = a*(1-e^2)/(1+e*cos(nu));
rhat = [cos(O)*cos(nu+w)-sin(O)*sin(nu+w)*cos(I);
        sin(O)*cos(nu+w)+sin(nu+w)*cos(I)*cos(O);
        sin(I)*sin(nu+w)];
dr_a = (1-e^2)/(1+e*cos(nu))*rhat*da;
dr_e = (-2*a*e*(1+e*cos(nu))-a*(1-e^2)*cos(nu))/(1+e*cos(nu)^2)*rhat*de;
dr_I = r*[sin(O)*sin(nu+w)*sin(I);
          -cos(O)*sin(nu+w)*sin(I);
          sin(nu+w)*cos(I)]*dI;
dr_O = r*[-sin(O)*cos(nu+w)-cos(O)*sin(nu+w)*cos(I);
          cos(O)*cos(nu+w)-sin(O)*sin(nu+w)*cos(I);
          0]*dO;
dr_w = r*[-cos(O)*sin(nu+w)-sin(O)*cos(nu+w)*cos(I);
          -sin(O)*sin(nu+w)+cos(nu+w)*cos(I)*cos(O);
          sin(I)*cos(nu+w)]*dw;
delta_r = dr_a+dr_e+dr_I+dr_O+dr_w;

end

