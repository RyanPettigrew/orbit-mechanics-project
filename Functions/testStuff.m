% Sat 1
%coe = [h, e, RA, inc, AoP, TA, a];
rEarth = 6378; % km
global mu
mu = 398600; % km^3/s^2

Sat1.e = 0; % Assume circular
Sat1.inc = deg2rad(0.8204);
Sat1.RAAN = deg2rad(83.0952);
Sat1.AoP = deg2rad(168.4676);
Sat1.epoch = 100; % Fix this
Sat1.TA = MA2TA(deg2rad(136.2128), Sat1.e);
Sat1.a = (35755 + rEarth + 35822 + rEarth)/2;
Sat1.h = sqrt(mu*(1 + Sat1.e)*Sat1.a); % Assume circular
[r,v] = coes2vector([Sat1.h, Sat1.e, Sat1.RAAN, Sat1.inc, Sat1.AoP, Sat1.TA, Sat1.a]);
Sat1.r = r;
Sat1.v = v;

plotOrbit(Sat1.r, Sat1.v, [0 24*60*60])

