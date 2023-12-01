
%% EXAMPLE: SCRIPT FOR CALLING FUNCTIONS/ORBIT CALULATIONS
clear; clc; close all;

run orbit_data.m
% run Time_orbits_project.m
run Orbits_propagations_orbits_project.m
run first_object_and_transfer.m
%run all_transfers.m


% These both do the same thing
disp(OBJECT1); 
disp(DATA.OBJECT1);

% How to get the data out of the structure; these both do the same thing
disp(OBJECT1.ECCENTRICITY); 
disp(DATA.OBJECT1.ECCENTRICITY);

