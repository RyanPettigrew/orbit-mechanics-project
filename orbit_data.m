%{
THIS FUNCTION LOADS THE .JSON ORBIT DATA STRUCTURE
EXAMPLE ON HOW TO USE IT IN ANY SCRIPT:
run orbit_data.m 
x = DATA.OBJECT1.ECCENTRICITY;
%}

clear; clc; close all;
myJsonFile = "orbits.json";
text = fileread(myJsonFile);
DATA = jsondecode(text);
