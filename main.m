
%% EXAMPLE: SCRIPT FOR CALLING FUNCTIONS/ORBIT CALULATIONS
clear; clc; close all;

run orbit_data.m

% These both do the same thing
disp(OBJECT1); 
disp(DATA.OBJECT1);

% How to get the data out of the structure; these both do the same thing
disp(OBJECT1.ECCENTRICITY); 
disp(DATA.OBJECT1.ECCENTRICITY);

