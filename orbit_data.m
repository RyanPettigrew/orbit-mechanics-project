%{
THIS FUNCTION LOADS THE .JSON ORBIT DATA STRUCTURE
EXAMPLE ON HOW TO USE IT IN ANY SCRIPT:
run orbit_data.m 
x = DATA.OBJECT1.ECCENTRICITY;
OR
run orbit_data.m 
x = OBJECT1.ECCENTRICITY;
%}

clear; clc; close all;
myJsonFile = "orbits.json";
text = fileread(myJsonFile);
DATA = jsondecode(text);
OBJECT1 = DATA.OBJECT1;
OBJECT2 = DATA.OBJECT2;
OBJECT3 = DATA.OBJECT3;
OBJECT4 = DATA.OBJECT4;
