%desattest.m
%Nathan Zimmerberg (nhz2)
%28SEP2018
%PAN

%Script to test time to desaturate reaction wheels using magnetorquer
%Run desatsim.slx simulink model with varius parameters and display 
%maxtime, the largest simulation end time, and i the index of that sim in simOutputs 
%Run the simulation again with the parameters that resulted in the
%largest end time.

%Constants
HMAG= 0.00675; % magnitude of initial angular momentum of the sat 
%with respect to sat center of mass, units kg*m^2/s
RMIN= 6771000; %minimum orbital radius, units m
RMAX= 6771000; %maximum orbital radius, units m
INCLMIN= 0.9; %minimum orbital inclination, units radians
INCLMAX= 0.9; %maximum orbital inclination, units radians
MOMENT= 0.08; %magnetic moment mag, units A*m^2
TIMEMAX= 24*3600;%inital starting time of the simulation, units seconds
N= 10; %number of simulations to run



%model parameters to vary
    %hCom_sat_eci initial angular momentum of the sat 
    %with respect to sat center of mass, units kg*m^2/s
    %vary on a sphere with radius HMAG
    
    %time_offset inital starting time of the simulation, units seconds
    %vary from 0s to 24*3600s incase earth rotates in a bad way
    
    %h_orbit_eci orbit angular momuntum/mass in ECI, units m^2/s
    %vary for a circular orbit from RMIN to RMAX radius
    
    %ri_orbit_eci intial orbital position in ECI,units m
    %vary for a circular orbit from RMIN to RMAX radius
    
    
    
%Script

[r, v]= orb2rv_s(RMIN,0.0,INCLMIN,0.0,0.0,0.0,3.986004418e14);
ri_orbit_eci= r;
h_orbit_eci = cross(r,v); %angular momuntum/mass in ECI units m^2/s

clear simIn
clear ri_orbit_ecis
clear h_orbit_ecis
clear time_offsets
clear hCom_sat_ecis
for i = 1:N
    rad= random('Uniform',RMIN,RMAX);
    incl= random('Uniform',INCLMIN,INCLMAX);
    [r, v]= orb2rv_s(rad,0.0,incl,0.0,0.0,0.0,3.986004418e14);
    ri_orbit_ecis(1:3,i)= r;
    h_orbit_ecis(1:3,i)= cross(r,v); %angular momuntum/mass in ECI units m^2/s7
    time_offsets(i)= random('Uniform',0.0,TIMEMAX);
    rand= random('Uniform',0.0,1,3,1);
    hCom_sat_ecis(1:3,i)= HMAG*rand/norm(rand);
    simIn(i)= Simulink.SimulationInput('desatsim');
    simIn(i)= setVariable(simIn(i),'hCom_sat_eci',hCom_sat_ecis(1:3,i));
    simIn(i)= setVariable(simIn(i),'time_offset',time_offsets(i));
    simIn(i)= setVariable(simIn(i),'h_orbit_eci',h_orbit_ecis(1:3,i));
    simIn(i)= setVariable(simIn(i),'ri_orbit_eci',ri_orbit_ecis(1:3,i));
end

simOutputs = sim(simIn);
% find the longest sim
clear time
for i = 1:N
    time(i)= simOutputs(i).getSimulationMetadata().ModelInfo.StopTime;
end
[maxtime,i]= max(time)

%rerun longest sim to get graphs
maxsim= Simulink.SimulationInput('desatsim');
maxsim= setVariable(simIn(i),'hCom_sat_eci',hCom_sat_ecis(1:3,i));
maxsim= setVariable(simIn(i),'time_offset',time_offsets(i));
maxsim= setVariable(simIn(i),'h_orbit_eci',h_orbit_ecis(1:3,i));
maxsim= setVariable(simIn(i),'ri_orbit_eci',ri_orbit_ecis(1:3,i));
ri_orbit_eci= ri_orbit_ecis(1:3,i);
h_orbit_eci= h_orbit_ecis(1:3,i);
time_offset= time_offsets(i);
hCom_sat_eci= hCom_sat_ecis(1:3,i);
maxsimOut = sim(maxsim);
