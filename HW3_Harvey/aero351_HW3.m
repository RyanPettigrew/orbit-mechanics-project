% Harvey Perkins


%% 4.15
clc
clear all

disp("4.15:")

mu = 398600; % km^3/s^2
rearth = 6378; % km

e = 1.5;
zp = 300; % km
rp = zp + rearth;
i = deg2rad(35); % rad
raan = deg2rad(130); % rad
aop = deg2rad(115); % rad
% At perigee, so
theta = 0; % rad

% get h
h = sqrt(mu*rp*(1 + e));

% r, v in perifocal frame
r = h^2/(mu*(1 + e))*[cos(theta);sin(theta);0];
disp("The radius vector in the perifocal frame is")
disp(r)
disp("km")
% As expected since theta = 0, r = rp

v = mu/h * [-sin(theta); e + cos(theta); 0];

disp("The velocity vector in the perifocal frame is")
disp(v)
disp("km/s")
% As expected, LEO alt but relatively eccentric

% Convert to ECI frame
CperitoECI = Cz(-raan)*Cx(-i)*Cz(-aop);

rECI = CperitoECI*r;

disp("The radius vector in the ECI frame is")
disp(rECI)
disp("km")

vECI = CperitoECI*v;

disp("The velocity vector in the ECI frame is")
disp(vECI)
disp("km/s")

% Basically impossible to heart check these


%% 5.6
clear all

disp("5.6:")

mu = 398600; % km^3/s^2
tol = 1e-8;

r1 = [5644; 2830; 4170]; % km
r2 = [-2240; 7320; 4980]; % km
dt = 20*60; % sec

% find deltatheta, assume prograde
cross = cross(r1,r2);
%if cross(3) >= 0
%    deltatheta = acos(dot(r1,r2)/(norm(r1)*norm(r2)));
%else
%    deltatheta = 2*pi - acos(dot(r1,r2)/(norm(r1)*norm(r2)));
%end

deltatheta = acos(dot(r1,r2)/(norm(r1)*norm(r2)));
deltatheta = asin(1*sqrt(1-cos(deltatheta)^2));

% get A
A = sqrt(norm(r1)*norm(r2))*sin(deltatheta)/sqrt(1 - cos(deltatheta));
if A == 0
    disp("Lambert's solution failed")
end

% Initial guess for z, bounds
z = 0;
zupper = 4*pi^2;
zlower = -zupper;
dt_calc = 0;

% precompute factorials
%k = 1:9;
%global factorial_C factorial_S
%factorial_C = (-1).^k./factorial(2*k + 2);
%factorial_S = (-1).^k./factorial(2*k + 3);

while abs(dt - dt_calc) > tol
    S = stumpfS_trig(z);
    C = stumpfC_trig(z);
    Y = norm(r1) + norm(r2) + A*(z*S - 1)/sqrt(C);
    chi = sqrt(Y/C);
    dt_calc = chi^3*S/sqrt(mu) + A*sqrt(Y)/sqrt(mu);
    if dt_calc < dt
        zlower = z;
    else
        zupper = z;
    end
    z = (zupper + zlower)/2;
end

% Calc final values
S = stumpfS_trig(z);
C = stumpfC_trig(z);
Y = norm(r1) + norm(r2) + A*(z*S - 1)/sqrt(C);
chi = sqrt(Y/C);
dt_calc = chi^3*S/sqrt(mu) + A*sqrt(Y)/sqrt(mu);

f = 1 - chi^2*C/norm(r1);
g = dt - chi^3*S/sqrt(mu);

fdot = sqrt(mu)*chi*(z*S - 1)/(norm(r1)*norm(r2));
gdot = 1 - chi^2*C/norm(r2);

v1 = (r2 - f*r1)/g;

r2_test = f*r1 + g*v1; % Heart check, works

v2 = fdot*r1 + gdot*v1;

disp("V1:")
disp(v1)
disp("V2:")
disp(v2)

%% 6.8
clear all

mu = 398600; % km^3/s^2
rearth = 6378; % km
% Circular to circular
r1 = 300 + rearth; % km
r3 = 3000 + rearth; % km

% velocities for orbits
v1 = sqrt(mu/r1);
v3 = sqrt(mu/r3);

% transfer eccentricity
e = (r3 - r1)/(r3 + r1);

% h transfer
h = sqrt(r1*mu*(1 + e*cos(0)));

% vp transfer
vp2 = h/r1;

% va transfer
va2 = h/r3;

% total dv
totaldv = vp2 - v1 + v3 - va2;
disp("6.8:");
disp("Total dv: " + totaldv + " km/s")

% Period of transfer orbit
a2 = (r1 + r3)/2;
T = 2*pi*a2^(3/2)/sqrt(mu);
disp("Transfer time: " + T/(2*60) + " min")

%% 6.23
clear all

mu = 398600; % km^3/s^2
rearth = 6378; % km

rp1 = 8100;
ra1 = 18900;
thetaB = deg2rad(45);
thetaC = deg2rad(150);

% Time for s/c C to reach point B
% Period - t150 + t45
a1 = (rp1 + ra1)/2; % SMA orbit 1
T1 = 2*pi*a1^(3/2)/sqrt(mu); % period orbit 1
e1 = (ra1 - rp1)/(ra1 + rp1);
% heart check: eccentric
h1 = sqrt(rp1*mu*(1 + e1));
t150 = TA2t(deg2rad(150),e1,mu,h1);
t45 = TA2t(deg2rad(45),e1,mu,h1);
% Heart check: 260 min reasonable for MEO/high LEO
% 150 degrees <-> 112 minutes less than half period
% 45 deg <-> 18 minutes pretty short

% Period of orbit 2
T2 = T1 - t150 + t45;

% Pe orbit 2
rp2 = h1^2/(mu*(1 + e1*cos(thetaB)));

% sma orbit 2
a2 = (T2*sqrt(mu)/(2*pi))^(2/3);

% Ap orbit 2
ra2 = 2*a2 - rp2;

% e orbit 2
e2 = (ra2 - rp2)/(ra2 + rp2);

% h orbit 2
h2 = sqrt(rp2*mu*(1 + e2));

% radial v orbit 1 @ B
vrad1B = mu*e1*sin(thetaB)/h1;

% perpendicular v orbit 1 @ B
vperp1B = mu*(1 + e1*cos(thetaB))/h1;

% perpendicular v orbit 2 @ Pe
vperp2 = h2/rp2; % wrong

% flight path angle orbit 1
gamma1 = atan2(e1*sin(thetaB),(1 + e1*cos(thetaB)));

% flight path angle orbit 2
gamma2 = 0;

% deltaV
dv = sqrt((vrad1B^2 + vperp1B^2) + vperp2^2 - 2*sqrt(vrad1B^2 + vperp1B^2)*vperp2*cos(gamma1 - gamma2));

dvTotal = 2*dv;
disp("6.23:")
disp("The total dv is: " + dvTotal + " km/s")


%% 6.25
clear all

mu = 398600; % km^3/s^2
rearth = 6378; % km
TA = deg2rad(100);

% orbit 1 stuff
rp1 = 1270 + rearth; % km
vp1 = 9; % km/s
h1 = rp1*vp1;
e1 = (h1^2/(mu*rp1) - 1);

% orbit 2 stuff
e2 = 0.4;
r2_100 = h1^2/(mu*(1 + e1*cos(TA)));
h2 = sqrt(r2_100*mu*(1 + e2*cos(TA)));

% flight path angle
gamma1 = atan2(e1*sin(TA), 1 + e1*cos(TA));
gamma2 = atan2(e2*sin(TA), 1 + e2*cos(TA));

deltagamma = gamma2 - gamma1;

% velocities
v1 = sqrt((mu*(1 + e1*cos(TA))/h1)^2 + (mu*e1*sin(TA)/h1)^2);
v2 = sqrt((mu*(1 + e2*cos(TA))/h2)^2 + (mu*e2*sin(TA)/h2)^2);

dv = sqrt(v1^2 + v2^2 - 2*v1*v2*cos(deltagamma));

disp("6.25:")
disp("The change in flight path angle is: " + rad2deg(deltagamma) + " deg")
disp("The total delta v is: " + dv + " km/s")

%% 6.44
clear all

mu = 398600; % km^3/s^2
rearth = 6378; % km

% Circular orbits
r1 = 300 + rearth;
r3 = 600 + rearth;
h1 = sqrt(r1*mu);
h2 = sqrt(r2*mu);
v1 = h1/r1;
v2 = h2/r2;

% a
e2 = (r3 - r1)/(r3 + r1);
h2 = sqrt(r1*mu*(1 + e2));
vp2 = h2/r1;
va2 = h2/r3;









function output = Cx(theta)
    output = [1,0,0;0,cos(theta),sin(theta);0,-sin(theta),cos(theta)];
end

function output = Cz(theta)
    output = [cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
end

% Stumpf function trig formulations
function C = stumpfC_trig(z)
    if z > 0
        C = (1 - cos(sqrt(z)))/z;
    elseif z < 0
        C = (cosh(sqrt(-z)) - 1)/(-z);
    else
        C = 1/2;
    end
end
function S = stumpfS_trig(z)
    if z > 0
        S = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
    elseif z < 0
        S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    else
        S = 1/6;
    end
end

function C = stumpfC(z)
    global factorial_C
    C = 1/2;
    for i = 1:9
        C = C + (z^i)*factorial_C(i);
    end
end
function S = stumpfS(z)
    global factorial_S
    S = 1/6;
    for i = 1:9
        S = S + (z^i)*factorial_S(i);
    end
end

function t = TA2t(TA,e,mu,h)
    % Elliptical orbits
    % TA in rad
    
    if (TA <= pi)
        E = 2*atan2(sqrt((1-e)/(1+e))*tan(TA/2),1);
    else
        E = 2*atan2(-sqrt((1-e)/(1+e))*tan(TA/2),-1);
    end
    
    Me = E - e*sin(E);
    
    T = 2*pi*(h/sqrt(1-e^2))^3/mu^2;
    
    t = Me*T/(2*pi);  
end