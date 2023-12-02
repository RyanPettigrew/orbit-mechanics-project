clc
clear all

global mu
mu = 398600; % km3/s2

rvect_object1 = [-2.77e4;2.26e4;0.0448e4];
vvect_object1 = [-2.11;-2.59;0.0265];
epoch = 8.856990254250000e+09;
% Propagate forwards 5 periods plus a bit
% period object 1 = 6.738e4 seconds
loiter_time_object1 = 5*5.738e4 + 2*60*60;
object1_depart_time = epoch + loiter_time_object1;
[rvect_object1_depart, vvect_object1_depart] = propagateOrbit(rvect_object1,vvect_object1,epoch,object1_depart_time);
%Plot Orbit #1
plotOrbit(rvect_object1, vvect_object1, [0 loiter_time_object1]);

% This Sim gives orbit #1
% simulateOrbit(rvect_object1, vvect_object1, [0 loiter_time_object1], rvect_object1_depart');

%% Lamberts to object 2
% At epoch
rvect_object2_new = [-0.0737e4; 3.10e4; 0.00863e4];
vvect_object2_new = [-3.62; -0.0609; -0.330];

lamberts12_time = 18*60*60;
object2_arrive_time = object1_depart_time + lamberts12_time;

[rvect_object2_arrive, vvect_object2_arrive] = propagateOrbit(rvect_object2_new,vvect_object2_new,epoch,object2_arrive_time);

% lamberts
[vsc_object1_depart, vsc_object2_arrive] = lamberts(rvect_object1_depart, rvect_object2_arrive, lamberts12_time, 1);

% Plot Orbit #2
plotOrbit(rvect_object1_depart, vsc_object1_depart, [0 lamberts12_time]);

% Propagate object 2 for 5 orbits
% period = 5.655e4
object2_wait_time = 5*5.655e4 + 30*60*60; % prev: 30
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
simulateOrbit(rvect_object2, vvect_object1, [0 loiter_time_object1], rvect_object1_depart');

%% Lamberts to object 3


lamberts23_time = 3.12*60*60; 
object3_arrive_time = object2_depart_time + lamberts23_time;

[rvect_object3_arrive, vvect_object3_arrive] = propagateOrbit(rvect_object3_new,vvect_object3_new,epoch,object3_arrive_time);

% lamberts
[vsc_object2_depart, vsc_object3_arrive] = lamberts(rvect_object2_depart, rvect_object3_arrive, lamberts23_time, 1);

plotOrbit(rvect_object2_depart, vsc_object2_depart, [0 lamberts23_time]);

% Propagate object 3 for 2 days
object3_depart_time = object3_arrive_time + 2*24*60*60;
[rvect_object3_depart, vvect_object3_depart] = propagateOrbit(rvect_object3_arrive,vvect_object3_arrive,object3_arrive_time,object3_depart_time);
plotOrbit(rvect_object3_arrive, vvect_object3_arrive, [0 2*24*60*60]);

disp("Object 2 departure burn: " + norm(vsc_object2_depart - vvect_object2_depart) + " km/s")
disp("Object 3 arrival burn: " + norm(vsc_object3_arrive - vvect_object3_arrive) + " km/s")


legend("Orbit 2", "Transfer 2-3", "Orbit 3")

hold off

%% Transfer to object 4

% Inc and RAAN change
coe_object3 = vector2coe(rvect_object3_depart, vvect_object3_depart,mu);
inc_object3 = coe_object3(4);
RAAN_object3 = coe_object3(3);
v_object3 = norm(vvect_object3_depart);

rvect_object4_start =  1.0e+02*[-4.299031324314775 1.923690022476296 7.777323084694023];
vvect_object4_start = [-17.907795963779567 -7.063735680315943 -8.198682991422665];
[rvect_object4, vvect_object4] = propagateOrbit(rvect_object4_start,vvect_object4_start,epoch,object3_depart_time);

coe_object4 = vector2coe(rvect_object4, vvect_object4,mu);
inc_object4 = coe_object4(4);
RAAN_object4 = coe_object4(3);

[delta_v] = inc_RAAN_change(v_object3,inc_object3,inc_object4,RAAN_object3,RAAN_object4);
% This calculates the new coes after the orbit at object 3 changes its inc
% and RAAN

coes_new = [coe_object3(1), coe_object3(2), coe_object4(3), coe_object4(4), coe_object3(5), coe_object3(6), coe_object3(7)];
[rvect_orbit3_inc_raan_change, vvect_orbit3_inc_raan_change] = coes2vector(coes_new,mu); % This is our new r and v vector that has an inc and RAAN the same as object 4 but everything else the same as object 3

% Hohmanns to get on same orbit
rp = norm(rvect_orbit3_inc_raan_change);
ra = norm(rvect_object4);
 [deltaVtotal,t] = hohmann(rp,ra);

time_after_hohmanns = object3_depart_time + t;

% HELP: how do I find the r and v vector? [rvect_posthohmanns,vvect_post_hohmanns] = r_and_v_of_hohmanns(rp,coes_new,coes_object4)
[rvect_object4_orbit,vvect_object4_orbit] = r_and_v_of_hohmanns(rvect_object3,rvect_object4);
[rvect_object4_posthohmann, vvect_object4_posthohmann] = propagateOrbit(rvect_object4,vvect_object4,epoch,time_after_hohmanns);
% Use period of the fourth object and time since periapse passage
% then use to find Me and use fourth's objects eccentricty.
%TA = MA2TA(Me,e)
% Phasing maneuver
% [deltaVtotal,t] = phasing_maneuver(rp,ra,TA)



%% Functions
function [rPrime, vPrime] = propagateOrbit(r, v, epoch, endTime)
    % Harvey Perkins
    % Propagates orbit from r,v vectors at epoch to endTime
    global mu

    t0 = 0;
    tf = endTime - epoch;

    options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Relative tolerance is the step in the function 

    [t,y] = ode45(@EOM, [t0 tf], [r;v], options);

    rPrime = y(end, 1:3)';
    vPrime = y(end, 4:6)';

end

function dy = EOM(t, y)
    global mu
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
    global mu
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
end

function [rvect,vvect] = r_and_v_of_hohmanns(rvect_object3,rvect_object4)
% First burn
rvect_object3 = [-2.77;2.26;0.0448];
rvect_object4 = [-4.299031324314775;1.923690022476296;7.777323084694023];
mu_earth = 98600;
r1 = norm(rvect_object3);
r2 = norm(rvect_object4);
at = (r1 + r2)/2;
vt1 = sqrt(mu_earth*((2/r1) - (1/at)));
h = r1*vt1;
true_anomaly_t = pi;
rt = r2;
rvect_t = rt*[cos(true_anomaly_t),sin(true_anomaly_t),0]
vvect_t = (h/rt)*[-sin(true_anomaly_t),cos(true_anomaly_t),0]
% Final burn
v2 = sqrt(mu_earth*((2/r2) - (1/at)));
% at arrival true_anomaly = 0;
true_anomaly_2 = 0;
rvect_2 = r2*[cos(true_anomaly_2),sin(true_anomaly_2),0]
vvect_2 = (h/r2)*[-sin(true_anomaly_2),cos(true_anomaly_2),0]
rvect = rvect_2
vvect = vvect_2 + vvect_t
end

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
function simulateOrbit(r, v, tspan, debris_positions)
    % r, v are column vectors for initial position & velocity
    % tspan is the timespan [t0 tf] for the orbit
    % debris_positions is a matrix of all positions of debris
    % Get EOM
    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    [t,y] = ode45(@EOM, tspan, [r;v], options);

    % Set up spacecraft STL
    stlData = stlread('SpaceX Crew Dragon - 3486142\files\SpaceXCrewDragon.stl');
    vertices = stlData.Points;
    faces = stlData.ConnectivityList;
    scale = 75; 
    vertices = vertices * scale;
    
    % Set up video sim
    video = VideoWriter('orbit_simulation.mp4', 'MPEG-4');
    open(video);
    
    % Earth texture
    earthTexture = 'earth.png';
    [X, Y, Z] = sphere(50);
    earth_radius = 6371;
    globe = surf(earth_radius*X, earth_radius*Y, earth_radius*Z, 'EdgeColor', 'none');
    set(globe, 'FaceColor', 'texturemap', 'CData', imread(earthTexture));
    hold on;
  
    % Plot initial orbit & debris
    orbitLine = plot3(y(:,1),y(:,2),y(:,3), 'b');
    debrisDots = plot3(debris_positions(:,1), debris_positions(:,2), debris_positions(:,3), 'ro');
    spacecraftDot = plot3(y(1,1), y(1,2), y(1,3), 'gs', 'MarkerSize', 10);
    xlabel('x, km');
    ylabel('y, km');
    zlabel('z, km');
    axis equal;
    grid on;

    % spacecraft
    spacecraftModel = patch('Vertices', vertices, 'Faces', faces, 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none');
    for spacecraft_idx = 1:length(t)
        currentPosition = y(spacecraft_idx, 1:3);
        transformedVertices = bsxfun(@plus, vertices, currentPosition - vertices(1,:));
        set(spacecraftModel, 'Vertices', transformedVertices);
        
        drawnow;frame = getframe(gcf); % video simulation
        writeVideo(video, frame);
        pause(0.05); 
    end
    close(video);
    hold off;
end

function simulateOrbitsALL(all_r, all_v, all_tspans, debris_positions)
    % r, v are column vectors for initial position & velocity
    % tspan is the timespan [t0 tf] for the orbit
    % debris_positions is a matrix of all positions of debris
    % Get EOM

    % Earth texture
    earthTexture = 'earth.png';
    [X, Y, Z] = sphere(50);
    earth_radius = 6371;
    globe = surf(earth_radius*X, earth_radius*Y, earth_radius*Z, 'EdgeColor', 'none');
    set(globe, 'FaceColor', 'texturemap', 'CData', imread(earthTexture));
    hold on;
    
    % Set up video sim
    video = VideoWriter('all_orbit_simulation.mp4', 'MPEG-4');
    open(video);

    % Set up spacecraft STL
    stlData = stlread('SpaceX Crew Dragon - 3486142\files\SpaceXCrewDragon.stl');
    vertices = stlData.Points;
    faces = stlData.ConnectivityList;
    scale = 25; 
    vertices = vertices * scale;
    spacecraftModel = patch('Vertices', vertices, 'Faces', faces, 'FaceColor', [0.8, 0.8, 0.8], 'EdgeColor', 'none');

    % All orbits and transfers
    for segment_idx = 1:length(all_r)
        r = all_r{segment_idx};
        v = all_v{segment_idx};
        tspan = all_tspans{segment_idx};
        options = odeset('RelTol',1e-8,'AbsTol',1e-8);
        [t,y] = ode45(@EOM, tspan, [r;v], options);

        plot3(y(:,1),y(:,2),y(:,3), 'b');
        
        % spacecraft positions
        for spacecraft_idx = 1:length(t)
            currentPosition = y(spacecraft_idx, 1:3);
            transformedVertices = bsxfun(@plus, vertices, currentPosition - vertices(1,:));
            set(spacecraftModel, 'Vertices', transformedVertices);

            drawnow;
            frame = getframe(gcf);
            writeVideo(video, frame);
            pause(0.05); %video speed
        end
    end

    % Plot debris positions
    plot3(debris_positions(:,1), debris_positions(:,2), debris_positions(:,3), 'ro');
    xlabel('x, km'); ylabel('y, km'); zlabel('z, km');
    axis equal; grid on; hold off;
    close(video);
end
