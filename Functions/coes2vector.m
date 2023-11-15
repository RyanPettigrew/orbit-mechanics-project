function [r, v] = coes2vector(coe)
    % Harvey Perkins
    % Converts COEs to r v vector
    %coe = [h, e, RA, inc, w, TA, a];
    global mu
    h = coe(1);
    e = coe(2);
    RAAN = coe(3);
    inc = coe(4);
    w = coe(5);
    TA = coe(6);
    a = coe(7);

    % Perifocal frame
    rPERI = h^2/(mu*(1 + e*cos(TA)))*[cos(TA);sin(TA);0];
    vPERI = mu/h*[-sin(TA);e + cos(TA); 0];

    CECItoPERI = Cz(w)*Cx(inc)*Cz(RAAN);
    CPERItoECI = CECItoPERI';

    r = CPERItoECI*rPERI;
    v = CPERItoECI*vPERI;

end