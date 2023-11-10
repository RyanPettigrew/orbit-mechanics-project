%% Time Code

% Inputs
% Y = ;
% M = ;
% D = ;

% hour = 4;
% min = 0;
% sec = 0;

% Auriel's Code for Julian Date: 
[JD] = julian_date(Y,M,D,hour,min,sec);
disp('The Julian Date is ' + string(JD) + ' days')

% Matlab Comparison for Julian Date:
time = ["2023-09-022 4:00:00"];
data = datetime(time);
JD_matlab_function = juliandate(data);
disp('The Julian Date using the built in Matlab function is ' + string(JD_matlab_function) + ' days')

% Local Sidereal Time:
% Y = 2007;
% M = 12;
% D = 21;
% hour = 10;
% min = 0;
% sec = 0;

% east_long_degrees = 144;
% east_long_minutes = 58;
% east_long_seconds = 0;
% east_longitude = east_long_degrees + east_long_minutes/60 + east_long_seconds/3600;

[LST] = local_sidereal(Y,M,D,hour,min,sec,east_longitude);
disp('The LST for part a of question 2 is ' + string(LST) + ' degrees')

%% General ODE to get position at different times

% rvect = [3207 5459 2714];
% vvect = [-6.532 0.7835 6.142];
% timespan = [0 24*60*60]; % CONVERT 24 hours to seconds but go back to hours at the end!

muearth = 398600; % km3/s2

% Set to call ODE
initialstate = [rvect vvect]; % state vector
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Relative tolerance is the step in the function 

% Output = ode45(@fun, timespan, state, options, anything else)
[tnew,statenew] = ode45(@airplane, timespan, initialstate, options, muearth);

r_new = norm([statenew(end,1) statenew(end,2) statenew(end,3)]); % New position magnitude
v_new = norm([statenew(end,4) statenew(end,5) statenew(end,6)]);% New velocity magnitude

disp('The position vector magnitude  after 24 hours is ' + string(r_new) + ' km')
disp('The velocity of the position vector after 24 hours is ' + string(v_new) + ' km/s')

figure
plot(statenew(:,1), statenew(:,2)) % 2D plot of the orbit's trajectory
xlabel('x,km')
ylabel('y,km')
title('2D Plot of the Orbits Trajectory')

hold on
plot(statenew(1,1),statenew(1,2),'.','MarkerSize',20)
plot(statenew(end,1),statenew(end,2),'o','MarkerSize',10)
legend('Plot of the trajectory','Position at t=0','Position at t=24hrs','Location','southeast')
hold off

figure
plot3(statenew(:,1), statenew(:,2),statenew(:,3)) % 3D plot of the orbit's trajectory
xlabel('x,km')
ylabel('y,km')
zlabel('z,km')
title('3D Plot of the Orbits Trajectory')

hold on
plot3(statenew(1,1),statenew(1,2),statenew(1,3),'.','MarkerSize',20)
plot3(statenew(end,1),statenew(end,2),statenew(end,3),'o','MarkerSize',10)
legend('Plot of the trajectory','Position at t=0','Position at t=24hrs','Location','southeast')
hold off

%% Transfer Times:

%%  General Period eq:
a = (rat + rpt)/2;
T = ((2*pi)/sqrt(mu_earth))*a^(3/2);
disp('The transfer time is ' + string(t))

%% Hohman: 
t = T/2; % Transfer time [s]

%% Bi elliptical:
% Get semimajor axis of transfer orbits 2 and 3
% Use the following eq to get time of flight for two semiellipses of
% bielliptic transfer
t_bielliptical = (1/2)*((2*pi/sqrt(mu_earth))*a_2^3/2 + (2*pi/sqrt(mu_earth))*a_3^(3/2));

%% Phasing maneuvers:

E_b = 2*atan(sqrt((1-ecc1)/(1+ecc1))*tan(true_anomaly_b/2));
% Use Kepler's equation
t_ab = (T1/2*pi)*(E_b - ecc1*sin(E_b));


%% Functions

% Julian Date
function [JD] = julian_date(Y,M,D,hour,min,sec)
J_0 = 367*Y - floor((7*(Y + floor((M + 9)/12)))/4) + floor((275*M)/9) + D + 1721013.5; % Calc for Julian day number at 0 h UT
UT = hour + min/60 + sec/3600; % Universal standard time
JD = J_0 + (UT/24); % Final Julian date calc
end

% Local Sidereal
function [LST] = local_sidereal(Y,M,D,hour,min,sec,east_longitude)
J_0 = 367*Y - floor((7*(Y + floor((M + 9)/12)))/4) + floor((275*M)/9) + D + 1721013.5;
J2000 = 2451545; % time from January 1, 2000 at noon
julian_century = 36525; % days
T_0 = (J_0 - J2000)/julian_century;
theta_G0 = 100.4606184 + 36000.77004*T_0 + 0.000387933*T_0^2 - 2.58*(10^(-8))*T_0^3; %  Greenwich sidereal time at 0 hr UT

while theta_G0 > 360
    theta_G0 = theta_G0 - 360;
end 
while theta_G0 < 0
    theta_G0 + 360;
end
end

% General ODE

function dstate = airplane(time,state,muearth)

% Equation of motion for a two body system

x = state(1);
y = state(2);
z = state(3); 
dx = state(4);
dy = state(5);
dz = state(6);

r = norm([x y z]);

ddx = -muearth*x/r^3;
ddy = -muearth*y/r^3;
ddz = -muearth*z/r^3;

dstate = [dx;dy;dz;ddx;ddy;ddz]; % must be a column vector, integrating this will give the positions and velocities

end 






