function coe = vector2coe(R, V, mu)
%{
% Computes the classical orbital elements (COEs)
% Pass in the state vector (R,V) and mu
%}

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
