clc
clear all
close all

%global mu
mu = 398600; % km3/s2

% rvect_object1 = [-2.77e4;2.26e4;0.0448e4];
% vvect_object1 = [-2.11;-2.59;0.0265];
 rvect_object1 = [-3.267467850897391e4; 2.666846629075791e4; 0.052740529923902e4];
 vvect_object1 = [-1.941435016388355; -2.382720492211831; 0.024428097819707];

epoch = 8.856990254250000e+09;
% Propagate forwards 5 periods plus a bit
% period object 1 = 8.617e4 seconds
loiter_time_object1 = 5*8.617e4 + 9*60*60;
object1_depart_time = epoch + loiter_time_object1;
[rvect_object1_depart, vvect_object1_depart] = propagateOrbit(rvect_object1,vvect_object1,epoch,object1_depart_time);
figure
plotOrbit(rvect_object1, vvect_object1, [0 loiter_time_object1]);
hold on

%% Lamberts to object 2
% At epoch
 rvect_object2_new = [-0.058320776102072e4; 3.719852979087268e4; 0.007613527772417e4];
 vvect_object2_new = [-3.302738702424707; -0.029338075546502; 0.300849920350975];
% rvect_object2_new = [-0.0737e4; 3.10e4; 0.00863e4];
% vvect_object2_new = [-3.62; -0.0609; -0.330];

% [rvect_object2_new, vvect_object2_new] = propagateOrbit(rvect_object2_new,vvect_object2_new,epoch,object1_depart_time);

lamberts12_time = 20.68*60*60; % prev: 20.68
object2_arrive_time = object1_depart_time + lamberts12_time;

[rvect_object2_arrive, vvect_object2_arrive] = propagateOrbit(rvect_object2_new,vvect_object2_new,epoch,object2_arrive_time);

% lamberts
[vsc_object1_depart, vsc_object2_arrive] = lamberts(rvect_object1_depart, rvect_object2_arrive, lamberts12_time, 1);

plotOrbit(rvect_object1_depart, vsc_object1_depart, [0 lamberts12_time]);

% Propagate object 2 for 5 orbits
% period = 7.436e4
object2_wait_time = 5*7.436e4 + 11.3*60*60; % prev: 11.3
object2_depart_time = object2_arrive_time + object2_wait_time;
[rvect_object2_depart, vvect_object2_depart] = propagateOrbit(rvect_object2_arrive,vvect_object2_arrive,object2_arrive_time,object2_depart_time);
plotOrbit(rvect_object2_arrive, vvect_object2_arrive, [0 object2_wait_time]);

legend("Orbit 1", "Transfer 1-2", "Orbit 2")

hold off

disp("Object 1 departure burn: " + norm(vsc_object1_depart - vvect_object1_depart) + " km/s")
disp("Object 2 arrival burn: " + norm(vsc_object2_arrive - vvect_object2_arrive) + " km/s")

% Transfer 2
figure
plotOrbit(rvect_object2_arrive, vvect_object2_arrive, [0 object2_wait_time]);
hold on
%% Lamberts to object 3


lamberts23_time = 3.9*60*60; 
object3_arrive_time = object2_depart_time + lamberts23_time;
 rvect_object3 = [0.099221289746647e3; -2.495960051500151e3; -6.829798227857773e3];
 vvect_object3 = [7.309683122916486; 1.173315870778262; -0.331602416389682];
% rvect_object3 = 10^2*[8.762180894975575;2.048550193447096;1.347181696876256];
% vvect_object4 = [-4.437140387214239;6.361898900101488;19.423106726508873];

[rvect_object3_arrive, vvect_object3_arrive] = propagateOrbit(rvect_object3,vvect_object3,epoch,object3_arrive_time);

% lamberts
[vsc_object2_depart, vsc_object3_arrive] = lamberts(rvect_object2_depart, rvect_object3_arrive, lamberts23_time, 1);

plotOrbit(rvect_object2_depart, vsc_object2_depart, [0 lamberts23_time]);
hold on
% Propagate object 3 for 5 periods
% period object 3 = 6.190e3
object3_wait_time = 5*6.190e3;
object3_depart_time = object3_arrive_time + object3_wait_time;
[rvect_object3_depart, vvect_object3_depart] = propagateOrbit(rvect_object3_arrive,vvect_object3_arrive,object3_arrive_time,object3_depart_time);
plotOrbit(rvect_object3_arrive, vvect_object3_arrive, [0 object3_wait_time]);

disp("Object 2 departure burn: " + norm(vsc_object2_depart - vvect_object2_depart) + " km/s")
disp("Object 3 arrival burn: " + norm(vsc_object3_arrive - vvect_object3_arrive) + " km/s")


legend("Orbit 2", "Transfer 2-3", "Orbit 3")

hold off

%% Transfers to object 4

% Inc and RAAN change
coe_object3 = vector2coe(rvect_object3_depart, vvect_object3_depart,mu);
inc_object3 = coe_object3(4);
RAAN_object3 = coe_object3(3);
v_object3 = norm(vvect_object3_depart);

rvect_object4_start = [-6.821074905741718e3;-0.609979433471615e3;2.479233696115208e3];
vvect_object4_start = [-2.083352744891964;-2.878043758469060;-6.490413688087917];

[rvect_object4_start, vvect_object4_start] = propagateOrbit(rvect_object4_start,vvect_object4_start,epoch,object3_depart_time);

coe_object4 = vector2coe(rvect_object4_start, vvect_object4_start,mu);
inc_object4 = coe_object4(4);
RAAN_object4 = coe_object4(3);

[delta_v_from_inc_RAAN_change] = inc_RAAN_change(v_object3,inc_object3,inc_object4,RAAN_object3,RAAN_object4);
disp('Object 3 inc and RAAN change = ' + string(delta_v_from_inc_RAAN_change) + ' km/s')
% This calculates the new coes after the orbit at object 3 changes its inc
% and RAAN

coes_new = [coe_object3(1), coe_object3(2), coe_object4(3), coe_object4(4), coe_object3(5), coe_object3(6), coe_object3(7)];
[rvect_orbit3_inc_raan_change, vvect_orbit3_inc_raan_change] = coes2vector(coes_new); % This is our new r and v vector that has an inc and RAAN the same as object 4 but everything else the same as object 3
figure
plotOrbit(rvect_object3_depart,vvect_object3_depart, [0 2*24*60*60]);
hold on
plotOrbit(rvect_orbit3_inc_raan_change,vvect_orbit3_inc_raan_change,[0 2*24*60*60]);
legend("Orbit 3", "Orbit 3 w/ Inc RAAN change")
title('Inc and RAAN Change')
grid on
hold off

% Hohmanns to get on same orbit
rp = norm(rvect_orbit3_inc_raan_change);
ra = norm(rvect_object4_start);
 [deltaVtotal,t] = hohmann(rp,ra);
disp('Hohmanns to Object 4 orbit = ' + string(deltaVtotal) + ' km/s')

time_after_hohmanns = object3_depart_time + t;

[rvect_object4_orbit,vvect_object4_orbit,rvect_t,vvect_t] = r_and_v_of_hohmanns(rvect_orbit3_inc_raan_change,rvect_object4_start);
[rvect_object4_posthohmann, vvect_object4_posthohmann] = propagateOrbit(rvect_object4_start,vvect_object4_start,epoch,time_after_hohmanns);

figure
plotOrbit(rvect_orbit3_inc_raan_change,vvect_orbit3_inc_raan_change,[0 2*24*60*60]);
hold on
plotOrbit(rvect_object4_start,vvect_object4_start,[0 2*24*60*60]);
hold on
plotOrbit(rvect_t,vvect_t,[0 t]);
legend("Orbit 3 w/ Inc RAAN change","Orbit 4","Hohmann Transfer")
title('Hohmann Transfer')
grid on
hold off

coes_new_object4 = vector2coe(rvect_object4_posthohmann',vvect_object4_posthohmann',mu);
% Phasing maneuver
[deltaVtotal,t] = phasing_maneuver(rvect_object4_orbit,rvect_object4_posthohmann,coes_new_object4(6),mu);
disp('Phasing maneuver to Object 4 = ' + string(deltaVtotal) + ' km/s')

function [rPrime, vPrime] = propagateOrbit(r, v, epoch, endTime)
    % Harvey Perkins
    % Propagates orbit from r,v vectors at epoch to endTime
    mu = 398600; % km3/s2
    t0 = 0;
    tf = endTime - epoch;

    options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Relative tolerance is the step in the function 

    [t,y] = ode45(@EOM, [t0 tf], [r;v], options);

    rPrime = y(end, 1:3)';
    vPrime = y(end, 4:6)';

end

function dy = EOM(t, y)
    mu = 398600;
    % Do stuff
    % y = [x; y; z; vx; vy; vz]

    r = norm([y(1), y(2), y(3)]);
    v = norm([y(4), y(5), y(6)]);
    dy(1,1) = y(4);
    dy(2,1) = y(5);
    dy(3,1) = y(6);
    dy(4,1) = -mu*y(1)/r^3;
    dy(5,1) = -mu*y(2)/r^3;
    dy(6,1) = -mu*y(3)/r^3;
end
function plotOrbit(r, v, tspan)
    % Harvey Perkins
    % plots orbit starting at given r,v vectors over given timespan [t0 tf]
    % r v column vectors

    options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Relative tolerance is the step in the function 

    [t,y] = ode45(@EOM, tspan, [r;v], options);

    %figure
    plot3(y(:,1),y(:,2),y(:,3))
    xlabel("x, km")
    xlabel("y, km")
    xlabel("y, km")
    axis equal
    hold on
    grid on;

end

function [v1, v2] = lamberts(r1, r2, dt, sign)

    % Implement lambert's method to calc velocities on transfer
    mu = 398600;
    tol = 1e-8;

    % find deltatheta

    crossprod = cross(r1,r2);
    z = crossprod(3);

    if sign > 0
        if z < 0
            deltatheta = 2*pi - acos(dot(r1,r2)/norm(r1)/norm(r2));
        else
            deltatheta = acos(dot(r1,r2)/norm(r1)/norm(r2));
        end
    else
        if z < 0
            deltatheta = acos(dot(r1,r2)/norm(r1)/norm(r2));
        else
            deltatheta = 2*pi - acos(dot(r1,r2)/norm(r1)/norm(r2));
        end
    end
    
    %deltatheta = acos(dot(r1,r2)/(norm(r1)*norm(r2)));
    %deltatheta = asin(1*sqrt(1-cos(deltatheta)^2));

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

function coe = vector2coe(R, V, mu)
% Computes the classical orbital elements (COEs)
% Pass in the state vector (R,V) and mu
r = norm(R);
v = norm(V);
vr = dot(R,V)/r;
H = cross(R,V);
h = norm(H);
% Uses Eq 4.7:
inc = acos(H(3)/h);
% Uses Eq 4.8:
N = cross([0 0 1],H);
n = norm(N);
% Uses Eq 4.9:
if n ~= 0
    RA = acos(N(1)/n);
    if N(2) < 0
        RA = 2*pi - RA;
    end
else
    RA = 0;
end
%Eq 4.10:
E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
e = norm(E);
% Small number below ecc to be zero
eps_min = 1.e-10;
%Eq 4.12 (for case e = 0):
if n ~= 0
    if e > eps_min
        w = acos(dot(N,E)/n/e);
        if E(3) < 0
            w = 2*pi - w;
        end
    else
        w = 0;
    end
else
    w = 0;
end
% Uses Eq 4.13a (for case e = 0):
if e > eps_min
    TA = acos(dot(E,R)/e/r);
    if vr < 0
        TA = 2*pi - TA;
    end
else
    cp = cross(N,R);
    if cp(3) >= 0
        TA = acos(dot(N,R)/n/r);
    else
        TA = 2*pi - acos(dot(N,R)/n/r);
    end
end
%Eq 4.62 (when a < 0 for a hyperbola):
a=h^2/mu/(1 - e^2);
% Return elements in COE vector
coe = [h, e, RA, inc, w, TA, a];
end

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

function [delta_v] = inc_RAAN_change(v_object3,inc_object3,inc_object4,RAAN_object3,RAAN_object4)
inc_initial = inc_object3;
inc_final = inc_object4;
RAAN_initial = RAAN_object3;
RAAN_final = RAAN_object4;

inc_change = abs(inc_initial - inc_final);
RAAN_change = abs(RAAN_initial - RAAN_final);

alpha = acos(cos(inc_initial)*cos(inc_final)+sin(inc_initial)*sin(inc_final)*cos(RAAN_change));
delta_v = 2*(v_object3)*sin(alpha/2);
end

function [deltaVtotal,t] = hohmann(rp,ra)
mu_earth = 398600; % km3/s2
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

deltaVi = abs(Vpt - v1);
deltaVf = abs(v3 - Vat);
deltaVtotal = deltaVi + deltaVf;

a = (rat + rpt)/2;
T = ((2*pi)/sqrt(mu_earth))*a^(3/2);
t = T/2; % Transfer time [s]
t = t/60; % Transfer time [min]
end

function [rvect,vvect,rvect_t,vvect_t] = r_and_v_of_hohmanns(rvect_object3,rvect_object4)
% First burn
mu_earth = 98600;
r1 = norm(rvect_object3);
r2 = norm(rvect_object4);
at = (r1 + r2)/2;
vt1 = sqrt(mu_earth*((2/r1) - (1/at)));
h = r1*vt1;
true_anomaly_t = pi;
rt = r2;
rvect_t = rt*[cos(true_anomaly_t),sin(true_anomaly_t),0];
vvect_t = (h/rt)*[-sin(true_anomaly_t),cos(true_anomaly_t),0];
% Final burn
v2 = sqrt(mu_earth*((2/r2) - (1/at)));
% at arrival true_anomaly = 0;
true_anomaly_2 = 0;
rvect_2 = r2*[cos(true_anomaly_2),sin(true_anomaly_2),0];
vvect_2 = (h/r2)*[-sin(true_anomaly_2),cos(true_anomaly_2),0];
rvect = rvect_2;
vvect = vvect_2 + vvect_t;
end

% Currently editing
% phasing maneuver - reviewing 
function [deltaVtotal,t] = phasing_maneuver(r_posthohmann,r_object4,TA)
rp = r_posthohmann;
ra = r_object4;

h_old = sqrt(2*mu_earth)*sqrt((rp*ra)/(rp+ra));
ecc_old = (ra-rp)/(ra+rp);
a_old = (1/2)*(rp+ra);
T = ((2*pi)/sqrt(mu_earth))*a_old^(3/2);
true_anomaly = deg2rad(TA);
E_c = 2*atan(sqrt((1-ecc_old)/(1+ecc_old))*tan(true_anomaly/2));
t_bc = (T/(2*pi))*(E_c-ecc_old*sin(E_c));
T2 = T - t_bc;
a_new = (((sqrt(mu_earth))*T2)/(2*pi))^(2/3);
rd = 2*a_new - rp;
h_new = sqrt(2*mu_earth)*sqrt((rp*rd)/(rp+rd));

vp_new_orbit = h_new/rp;
vp_old_orbit = h_old/rp;

deltaStart = vp_new_orbit + vp_old_orbit;
deltaEnd = vp_old_orbit - vp_new_orbit;
deltaVtotal = abs(deltaStart) + abs(deltaEnd);

t = t_bc/60; % mins
end
