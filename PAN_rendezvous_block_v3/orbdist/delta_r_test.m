clear
close all

muE = 3.98600441500E+14;
nsteps = 10001;
nuspan = linspace(0,2*pi,nsteps);

termination_states

%% ISS rendezvous

[a_t_ISS, e_t_ISS, I_t_ISS, Om_t_ISS, w_t_ISS] = x2orb(x_t_ISS,muE);
els_t_ISS = [a_t_ISS; e_t_ISS; I_t_ISS; Om_t_ISS; w_t_ISS];
[a_f_ISS, e_f_ISS, I_f_ISS, Om_f_ISS, w_f_ISS] = x2orb(x_f_ISS,muE);
els_f_ISS = [a_f_ISS; e_f_ISS; I_f_ISS; Om_f_ISS; w_f_ISS];
d_els_ISS = els_t_ISS - els_f_ISS;
deltars_ISS = zeros(3,nsteps);
deltarns_ISS = zeros(1,nsteps);

for i=1:nsteps
    deltars_ISS(:,i) = delta_r(els_t_ISS,d_els_ISS,nuspan(i));
    deltarns_ISS(i) = norm(deltars_ISS(:,i));
end

%% 400-500km Hohmann entry
[a_t_h451, e_t_h451, I_t_h451, Om_t_h451, w_t_h451] = x2orb(x_t_h451,muE);
els_t_h451 = [a_t_h451, e_t_h451, I_t_h451, Om_t_h451, w_t_h451];
[a_f_h451, e_f_h451, I_f_h451, Om_f_h451, w_f_h451] = x2orb(x_f_h451,muE);
els_f_h451 = [a_f_h451, e_f_h451, I_f_h451, Om_f_h451, w_f_h451];
d_els_h451 = els_t_h451 - els_f_h451;
deltars_h451 = zeros(3,nsteps);
deltarns_h451 = zeros(1,nsteps);

for i=1:nsteps
    deltars_h451(:,i) = delta_r(els_t_h451,d_els_h451,nuspan(i));
    deltarns_h451(i) = norm(deltars_h451(:,i));
end

%% 400-500km Hohmann exit
[a_t_h452, e_t_h452, I_t_h452, Om_t_h452, w_t_h452] = x2orb(x_t_h452,muE);
els_t_h452 = [a_t_h452, e_t_h452, I_t_h452, Om_t_h452, w_t_h452];
[a_f_h452, e_f_h452, I_f_h452, Om_f_h452, w_f_h452] = x2orb(x_f_h452,muE);
els_f_h452 = [a_f_h452, e_f_h452, I_f_h452, Om_f_h452, w_f_h452];
d_els_h452 = els_t_h452 - els_f_h452;
deltars_h452 = zeros(3,nsteps);
deltarns_h452 = zeros(1,nsteps);

for i=1:nsteps
    deltars_h452(:,i) = delta_r(els_t_h452,d_els_h452,nuspan(i));
    deltarns_h452(i) = norm(deltars_h452(:,i));
end

%% ISS rendezvous 2

[a_t_ISS2, e_t_ISS2, I_t_ISS2, Om_t_ISS2, w_t_ISS2] = x2orb(x_t_ISS2,muE);
els_t_ISS2 = [a_t_ISS; e_t_ISS2; I_t_ISS2; Om_t_ISS2; w_t_ISS2];
[a_f_ISS2, e_f_ISS2, I_f_ISS2, Om_f_ISS2, w_f_ISS2] = x2orb(x_f_ISS2,muE);
els_f_ISS2 = [a_f_ISS2; e_f_ISS2; I_f_ISS2; Om_f_ISS2; w_f_ISS2];
d_els_ISS2 = els_t_ISS2 - els_f_ISS2;
deltars_ISS2 = zeros(3,nsteps);
deltarns_ISS2 = zeros(1,nsteps);

for i=1:nsteps
    deltars_ISS2(:,i) = delta_r(els_t_ISS2,d_els_ISS,nuspan(i));
    deltarns_ISS2(i) = norm(deltars_ISS2(:,i));
end

%% Test 1 - [405km, .001, 50.1, 189.9, 283] to [420km, .0001, 50, 190, 282]

[deltars_test1,deltarns_test1] = delta_r_data(x_t_test1,x_f_test1,muE,nsteps);

%% Plotting
figure(1)
plot(nuspan,deltarns_ISS)
hold on
plot(nuspan,deltarns_ISS2)
plot(nuspan,deltarns_test1)
plot(nuspan,deltarns_h451)
plot(nuspan,deltarns_h452)