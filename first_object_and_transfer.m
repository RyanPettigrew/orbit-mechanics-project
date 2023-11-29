clc
clear all

global mu
mu = 398600; % km3/s2

rvect_object1 = [-2.77e4; 2.26e4; 0.0448e4];
vvect_object1 = [-2.11; -2.59; 0.0265];
epoch = 8.856990254250000e+09;

% Propagate forwards 2 days
object1_depart_time = epoch + 2*24*60*60;
[rvect_object1_depart, vvect_object1_depart] = propagateOrbit(rvect_object1,vvect_object1,epoch,object1_depart_time);
plotOrbit(rvect_object1, vvect_object1, [0 2*24*60*60]);

% Lamberts to object 2
% At epoch
rvect_object2_new = [-0.0737e4; 3.10e4; 0.00863e4];
vvect_object2_new = [-3.62; -0.0609; -0.330];

lamberts12_time = 6*60*60;
object2_arrive_time = object1_depart_time + lamberts12_time;

[rvect_object2_arrive, vvect_object2_arrive] = propagateOrbit(rvect_object2_new,vvect_object2_new,epoch,object2_arrive_time);

% lamberts
[vsc_object1_depart, vsc_object2_arrive] = lamberts(rvect_object1_depart, rvect_object2_arrive, lamberts12_time, -1);
%[vsc_object1_depart, vsc_object2_arrive] = laberts_MW_that_works(rvect_object1_depart, rvect_object2_arrive, lamberts12_time, -1, mu);

plotOrbit(rvect_object1_depart, vsc_object1_depart, [0 lamberts12_time]);

% Propagate object 2 for 2 days
object2_depart_time = object2_arrive_time + 2*24*60*60;
[rvect_object2_depart, vvect_object2_depart] = propagateOrbit(rvect_object2_arrive,vvect_object2_arrive,object2_arrive_time,object2_depart_time);
plotOrbit(rvect_object2_arrive, vvect_object2_arrive, [0 2*24*60*60]);

legend("Orbit 1", "Transfer 1-2", "Orbit 2")

hold off



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

function [v_1,v_2] = laberts_MW_that_works(r_1,r_2,tspan,sign,mu)

% shall utilize lamberts method to find the velocity 

tol = 10^(-7) ;

% finding z of r1 x r2
rcross = cross(r_1,r_2);
z_12 = rcross(3);

% Calculate magnitudes of r1 and r2
r_1_mag = norm(r_1) ;
r_2_mag = norm(r_2) ;

% obtaining theta

if sign > 0
    if z_12 < 0
        theta = 360 - acosd(dot(r_1,r_2)/r_1_mag/r_2_mag);
    else
        theta = acosd(dot(r_1,r_2)/r_1_mag/r_2_mag);
    end
else
    if z_12 < 0
        theta = acosd(dot(r_1,r_2)/r_1_mag/r_2_mag);
    else
        theta = 360 - acosd(dot(r_1,r_2)/r_1_mag/r_2_mag);
    end
end

 A = sind(theta) * sqrt((r_1_mag*r_2_mag) / (1-cosd(theta))) ;

y_fun =@(z) r_1_mag + r_2_mag + A*((z*stumpS(z)-1)/(sqrt(stumpC(z)))) ;
w_fun =@(z) sqrt(y_fun(z)/stumpC(z)) ;
looper_fun =@(z) (((w_fun(z))^(3))*stumpS(z))/(sqrt(mu)) + (A*sqrt(y_fun(z)))/(sqrt(mu)) ;

% Bi-section

z_upper = 4*pi^2 ;
z_lower = -4*pi^2 ;
z = (z_upper+z_lower) /2 ;

y = y_fun(z) ;
w = w_fun(z) ;
looper = looper_fun(z) ;

while abs(looper - tspan) > tol
    z = (z_upper+z_lower) /2 ;
    y = y_fun(z) ;
    w = w_fun(z) ;
    looper = looper_fun(z) ;

    if looper <= tspan
        z_lower = z ;
    else
        z_upper = z;
    end
    
end

% happy face

f = 1 - y/r_1_mag ;
d_f = sqrt(mu*y/stumpC(z))*(z*stumpS(z)-1)/r_1_mag/r_2_mag ;
g = A*sqrt(y/mu) ;
d_g = 1 - y/r_2_mag ;

v_1 = (1/g)*(r_2 - (f*r_1)) ; 
v_2 = d_f*r_1+(d_g*v_1) ; 

end

% Stumpf function trig formulations
function C = stumpC(z)
    if z > 0
        C = (1 - cos(sqrt(z)))/z;
    elseif z < 0
        C = (cosh(sqrt(-z)) - 1)/(-z);
    else
        C = 1/2;
    end
end
function S = stumpS(z)
    if z > 0
        S = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
    elseif z < 0
        S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    else
        S = 1/6;
    end
end