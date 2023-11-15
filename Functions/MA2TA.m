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