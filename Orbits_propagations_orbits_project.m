%% Orbits project propagated Objects
% This code first gets the r and v vector from the the COEs and TLEs and
% then propagates them backwards to the epoch of the first object using ODE
% 45.
format long
clc; clear all; close all;
mu_earth = 398600; % km3/s2
global mu
mu = mu_earth;

%% Object 1: Start time is the epoch
% Julien Date of Object 1
% 26 November 2023 13:41:42
Y = 2023;
M = 11;
D = 26;
hour = 13;
min = 41;
sec = 42;

[JD] = julian_date(Y,M,D,hour,min,sec);
JD_seconds_object1 = JD*24*3600;
epoch = JD_seconds_object1;

% r and v vector to get COEs
%coe = [h, e, RA, inc, w, TA, a];
RA = deg2rad(83.2529);
e = 0.0007950;
inc = deg2rad(0.8492);
w = deg2rad(170.6011);
Me = deg2rad(247.0120);
TA = MA2TA(Me,e);
rp = 35755; % [km]
ra = 35822; % [km]
a = (ra + rp)/2; % Semi major axis
h = sqrt(a*mu_earth*(1-e^2));
coe = [h, e, RA, inc, w, TA, a];
[rvect_object1, vvect_object1] = coes2vector(coe);
r_object1 = norm(rvect_object1);
v_object1 = norm(vvect_object1);

%% Object 2
% Julien date of object 2
Y = 2023;
M = 11;
D = 26;
hour = 11;
min = 45;
sec = 52;
[JD] = julian_date(Y,M,D,hour,min,sec);
JD_seconds_object2 = JD*24*3600;
time_difference1 = JD_seconds_object2 - JD_seconds_object1;
% r and v vector to get COEs
%coe = [h, e, RA, inc, w, TA, a];
% epoch = 26 November 2023 11:45:52
RA = deg2rad(89.6110);
e = 0.0275217;
inc = deg2rad(5.2052);
w = deg2rad(346.2632);
Me = deg2rad(12.8256);
TA = MA2TA(Me,e);
rp = 30791; % [km]
ra = 32895; % [km]
a = (ra + rp)/2; % Semi major axis
h = sqrt(a*mu_earth*(1-e^2));
coe = [h, e, RA, inc, w, TA, a];
[rvect_object2, vvect_object2] = coes2vector(coe);

% Propagated backwards in time
% Output = ode45(@fun, timespan, state, options, anything else)
timespan = [time_difference1 0]; % CONVERT 24 hours to seconds but go back to hours at the end!
% Set to call ODE
rvect = rvect_object2;
vvect = vvect_object2;
initialstate = [rvect vvect]; % state vector
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Relative tolerance is the step in the function 
[tnew,statenew] = ode45(@airplane, timespan, initialstate, options, mu_earth);
rvect_object2_new = [statenew(end,1); statenew(end,2); statenew(end,3)];
vvect_object2_new = [statenew(end,4); statenew(end,5); statenew(end,6)];
r_object2 = norm([statenew(end,1) statenew(end,2) statenew(end,3)]); % New position magnitude
v_object2 = norm([statenew(end,4) statenew(end,5) statenew(end,6)]);% New velocity magnitude

%% Object 3: EGRS 3

% Julien date of object 3
Y = 2023;
M = 11;
D = 27;
hour = 13;
min = 51;
sec = 52;
[JD] = julian_date(Y,M,D,hour,min,sec);
JD_seconds_object3 = JD*24*3600;
time_difference2 = JD_seconds_object3 - JD_seconds_object1;
% r and v vector to get COEs
%coe = [h, e, RA, inc, w, TA, a];
% epoch = 27 November 2023 13:51:52
RA = deg2rad(10.0491);
e = 0.0022628;
inc = deg2rad(70.0805);
w = deg2rad(236.9837);
Me = deg2rad(241.0715);
TA = MA2TA(Me,e);
rp = 892; % [km]
ra = 925; % [km]
a = (ra + rp)/2; % Semi major axis
h = sqrt(a*mu_earth*(1-e^2));
coe = [h, e, RA, inc, w, TA, a];
[rvect_object3, vvect_object3] = coes2vector(coe);

% Propagated backwards in time
% Output = ode45(@fun, timespan, state, options, anything else)
timespan = [time_difference2 0]; % CONVERT 24 hours to seconds but go back to hours at the end!
% Set to call ODE
rvect = rvect_object3;
vvect = vvect_object3;
initialstate = [rvect vvect]; % state vector
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Relative tolerance is the step in the function 
[tnew,statenew] = ode45(@airplane, timespan, initialstate, options, mu_earth);
rvect_object3_new = [statenew(end,1); statenew(end,2); statenew(end,3)];
vvect_object3_new = [statenew(end,4); statenew(end,5); statenew(end,6)];
r_object3 = norm([statenew(end,1) statenew(end,2) statenew(end,3)]); % New position magnitude
v_object3 = norm([statenew(end,4) statenew(end,5) statenew(end,6)]);% New velocity magnitude

%% Object 4: GGSE 2 - Orbit

% Julien date of object 4
Y = 2023;
M = 11;
D = 27;
hour = 12;
min = 42;
sec = 53;
[JD] = julian_date(Y,M,D,hour,min,sec);
JD_seconds_object4 = JD*24*3600;
time_difference4 = JD_seconds_object4 - JD_seconds_object1;
% r and v vector to get COEs
%coe = [h, e, RA, inc, w, TA, a];
% epoch = 27 November 2023 12:42:53
RA = deg2rad(12.6485);
e = 0.0023500;
inc = deg2rad(70.0804);
w = deg2rad(239.6652);
Me = deg2rad(120.2130);
TA = MA2TA(Me,e);
rp = 891; % [km]
ra = 925; % [km]
a = (ra + rp)/2; % Semi major axis
h = sqrt(a*mu_earth*(1-e^2));
coe = [h, e, RA, inc, w, TA, a];
[rvect_object4, vvect_object4] = coes2vector(coe);

% Propagated backwards in time
% Output = ode45(@fun, timespan, state, options, anything else)
timespan = [time_difference4 0]; % CONVERT 24 hours to seconds but go back to hours at the end!
% Set to call ODE
rvect = rvect_object4;
vvect = vvect_object4;
initialstate = [rvect vvect]; % state vector
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Relative tolerance is the step in the function 
[tnew,statenew] = ode45(@airplane, timespan, initialstate, options, mu_earth);
rvect_object4_new = [statenew(end,1); statenew(end,2); statenew(end,3)];
vvect_object4_new = [statenew(end,4); statenew(end,5); statenew(end,6)];
r_object4 = norm([statenew(end,1) statenew(end,2) statenew(end,3)]); % New position magnitude
v_object4 = norm([statenew(end,4) statenew(end,5) statenew(end,6)]); % New velocity magnitude

%% Functions

% Julian Date
function [JD] = julian_date(Y,M,D,hour,min,sec)
J_0 = 367*Y - floor((7*(Y + floor((M + 9)/12)))/4) + floor((275*M)/9) + D + 1721013.5; % Calc for Julian day number at 0 h UT
UT = hour + min/60 + sec/3600; % Universal standard time
JD = J_0 + (UT/24); % Final Julian date calc
end

% Mean anomaly to true anomaly
function TA = MA2TA(Me,e)
    % Elliptical orbits
    
    tol = 1e-8;
    
    if (Me < pi)
        Eguess = Me + e/2;
    else
        Eguess = Me - e/2;
    end
    
    E(1) = Eguess;
    err = 10;
    
    while err > tol
        % Newton's method
        E(length(E) + 1) = E(end) - (Me - E(end) + e*sin(E))/(-1 + e*cos(E));
        err = abs(E(end) - E(end - 1));
    end
    
    E = E(end);
    
    % Unsure if this works for any angle
    if (Me <= pi)
        TA = 2*atan2(tan(E/2)/sqrt((1-e)/(1+e)),1); 
    else
        TA = 2*atan2(-tan(E/2)/sqrt((1-e)/(1+e)),-1);
    end

end

% COEs to vector
function [r, v] = coes2vector(coe)
    % Harvey Perkins
    % Converts COEs to r v vector
    %coe = [h, e, RA, inc, w, TA, a];
    mu_earth = 398600; % km3/s2
    h = coe(1);
    e = coe(2);
    RAAN = coe(3);
    inc = coe(4);
    w = coe(5);
    TA = coe(6);
    a = coe(7);

    % Perifocal frame
    rPERI = h^2/(mu_earth*(1 + e*cos(TA)))*[cos(TA);sin(TA);0];
    vPERI = mu_earth/h*[-sin(TA);e + cos(TA); 0];
    
    Cx_inc = [1,0,0;0,cos(inc),sin(inc);0,-sin(inc),cos(inc)];
    Cz_w = [cos(w),sin(w),0;-sin(w),cos(w),0;0,0,1];
    Cz_RAAN = [cos(RAAN),sin(RAAN),0;-sin(RAAN),cos(RAAN),0;0,0,1];
    CECItoPERI = Cz_w*Cx_inc*Cz_RAAN;
    CPERItoECI = CECItoPERI';

    r = CPERItoECI*rPERI;
    v = CPERItoECI*vPERI;

end

% ODE
function dstate = airplane(time,state,mu_earth)
% Equation of motion for a two body system
x = state(1);
y = state(2);
z = state(3); 
dx = state(4);
dy = state(5);
dz = state(6);

r = norm([x y z]);

ddx = -mu_earth*x/r^3;
ddy = -mu_earth*y/r^3;
ddz = -mu_earth*z/r^3;

dstate = [dx;dy;dz;ddx;ddy;ddz]; % must be a column vector, integrating this will give the positions and velocities
end 



