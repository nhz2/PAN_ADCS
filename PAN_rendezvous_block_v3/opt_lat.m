function lat_opt = opt_lat(r,v,x_t,x,K,L,muE)
%#codegen

z0=[r;v];

E=norm(v)^2/2-muE/norm(r);
a=-muE/(2*E);
T=2*pi/sqrt(muE)*a^(3/2);

options=odeset('AbsTol',1e-12,'RelTol',1e-12);
tspan=linspace(0,T,10000);
[~,zs] = ode45(@orbit1diff,tspan,z0,options,muE,[0;0;0]);

max_rat=0;
rat=0;
lat_opt=0;
N=length(tspan);
x_f=x;
x_diff=x_t-x_f;
% x_diff_hat=x_diff/norm(x_diff);

for i=1:N
    r_f=zs(i,1:3);
    v_f=zs(i,4:6);
    
%     if (abs(dot(r_f,v_f))<perpTol)
%         rat=-999999;
%     end
    
    A=[   0    -r_f(3)  r_f(2);
        r_f(3)    0    -r_f(1);
       -r_f(2)  r_f(1)    0   ;
        v_f(1)  v_f(2)  v_f(3)];

    B=-crs(v_f)*crs(r_f)/muE;    
    M=[A;
       B];
   
    rhat_dot=v_f'/norm(r_f)-dot(v_f,r_f)*r_f'/norm(r_f)^3;

    dV_vec=pinv(M)*(x_diff-[0;0;0;0;rhat_dot]);
%     dV_vec=pinv(M)*x_diff;
    dV_opt_sqd=dV_vec'*L*dV_vec;
    dx_vec=M*dV_vec+[0;0;0;0;rhat_dot];
    dx_opt_sqd=dx_vec'*K*dx_vec;
    
    rat=dx_opt_sqd/dV_opt_sqd;
    
    if (rat>max_rat)
        max_rat=rat;
        lat_opt=atan2(r_f(2),r_f(1));
        if (lat_opt<0)
            lat_opt=lat_opt+2*pi;
        end
    end
end
end
