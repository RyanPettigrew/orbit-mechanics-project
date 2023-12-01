function simulateOrbit(r_inital, v_inital, tspan, debris_positions)
    % r, v are column vectors for initial position & velocity
    % tspan is the timespan [t0 tf] for the orbit
    % debris_positions is a matrix where each row is a position vector of debris

    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    [t,y] = ode45(@EOM, tspan, [r_inital;v_inital], options);

    for spacecraft_idx = 1:length(t)
        clf; 
        plot3(y(:,1),y(:,2),y(:,3), 'b');
        hold on;

        % Spacecraft at current position
        spacecraft_position = y(spacecraft_idx, 1:3);
        drawSpacecraft(spacecraft_position);
        for i = 1:size(debris_positions, 1)
            plot3(debris_positions(i,1), debris_positions(i,2), debris_positions(i,3), 'ro')
        end
        % Draw Earth
        [X, Y, Z] = sphere;
        Earth_radius = 6371;
        surf(Earth_radius*X, Earth_radius*Y, Earth_radius*Z, 'EdgeColor', 'none');
        xlabel('x, km');
        ylabel('y, km');
        zlabel('z, km');
        axis equal;
        grid on;
        drawnow;
        pause(0.05);
    end
    hold off;
end
