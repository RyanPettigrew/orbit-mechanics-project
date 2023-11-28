%% Inc and RAAN Change
% If given inc and RAAN in radians
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