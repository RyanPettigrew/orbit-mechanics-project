clc
clear all

global mu
mu = 398600; % km3/s2

rvect_object1 = [-0.507e4; -3.496e4; 0.579e4];
vvect_object1 = [2.15; -0.72; -2.44];
epoch = 8.856990254250000e+09;

% Propagate forwards 2 days
dt = 2*24*60*60;
[rvect_object1_end, vvect_object1_end] = propagateOrbit(rvect_object1,vvect_object1,epoch,epoch + dt);
plotOrbit(rvect_object1, vvect_object1, [0 dt]);

% Object 2
rvect_object2_new = [-1.35e4; 1.59e4; -2.29e4];
vvect_object2_new = [-0.715; -3.095; -1.762];

% lamberts



function [rPrime, vPrime] = propagateOrbit(r, v, epoch, endTime)
    % Harvey Perkins
    % Propagates orbit from r,v vectors at epoch to endTime
    global mu

    t0 = 0;
    tf = endTime - epoch;

    options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Relative tolerance is the step in the function 

    [t,y] = ode45(@EOM, [t0 tf], [r;v], options);

    rPrime = y(end, 1:3);
    vPrime = y(end, 4:6);

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

    figure
    plot3(y(:,1),y(:,2),y(:,3))
    xlabel("x, km")
    xlabel("y, km")
    xlabel("y, km")

end
