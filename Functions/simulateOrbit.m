function simulateOrbit(r, v, tspan, debris_positions)
    % NO SPACECRAFT
    % r, v are column vectors for initial position & velocity
    % tspan is the timespan [t0 tf] for the orbit
    % debris_positions is a matrix where each row is a position vector of debris

    earthTexture = 'earth.png';
    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    [t,y] = ode45(@EOM, tspan, [r;v], options);

    %nEarth texture
    [X, Y, Z] = sphere(50); 
    radius_earth = 6371;
    globe = surf(radius_earth*X, radius_earth*Y, radius_earth*Z, 'EdgeColor', 'none');
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
    for spacecraft_idx = 1:length(t)
        set(spacecraftDot, 'XData', y(spacecraft_idx, 1), 'YData', y(spacecraft_idx, 2), 'ZData', y(spacecraft_idx, 3));
        drawnow;
        pause(0.05); 
    end

    hold off;
end
