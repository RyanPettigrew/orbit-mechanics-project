% rp = radius of perigee of original orbit
% ra = radius of apogee of new orbit

function [deltaVtotal,t] = hohmann(rp,ra)

rp = zp + r_earth; % Radius of perigee [km]
v1 = sqrt(mu_earth/rp); % Speed of og orbit; At LEO
h1 = rp*v1;

v3 = sqrt(mu_earth/ra); % Speed of final orbit; Should be slower than v1
h3 = ra*v3;

ecc = (ra - rp)/(rp + ra); % Should be an ellipse; eccentricity of the transfer orbit
rat = ra;
rf = ra;
rpt = rp;
h2 = sqrt(mu_earth*rpt*(1+ecc));
Vpt = h2/rpt; % Speed at perigee of transfer orbit
Vat = h2/rat; % Speed at apogee of transfer orbit; Slower than Vpt;

deltaVi = Vpt - v1;
deltaVf = v3 - Vat;
deltaVtotal = deltaVi + deltaVf;
disp('The total delta v is ' + string(deltaVtotal) + ' km/s')

a = (rat + rpt)/2;
T = ((2*pi)/sqrt(mu_earth))*a^(3/2);
t = T/2; % Transfer time [s]
t = t/60; % Transfer time [min]
disp('The transfer time is ' + string(t) + ' minutes')