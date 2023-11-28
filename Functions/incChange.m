% r1 = radius of original orbit
% r2 = radius of new orbit

function [tdelta,t] = incChange(r1,r2,deltainc)

mu = 398600;
v1 = sqrt(mu/r1);
v1_new = sqrt(mu/r2);
h = sqrt((2*mu*r1*r2)/(r1+r2));
v2 = h/r1;
v2_new = h/r2;

% calculate delta v's with a speed change and deltainc @ same time.
deltav1 = v2-v1;
deltav2 = sqrt(v2_new^2+v1_new^2-2*v2_new*v1_new*cosd(deltainc));
tdelta = deltav1*deltav2;

disp('The total delta v is ' + string(tdelta) + ' km/s')

a = (r1+r2)/2;
T = ((2*pi)/sqrt(mu))*a^(3/2);
t = (T/2)/60; %mins
disp('The transfer time is ' + string(t) + ' mins.')

end