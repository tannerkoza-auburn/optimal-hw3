%% Optimal Estimation - Homework 3 - Problem 3

clear 
clc
close all

%% Part A

data = readtable('hw3_3.txt');
data.Properties.VariableNames = {'t','E','N','psi','gyro','radar'};

figure
plot(data.t,data.E,'.')
title('East (m)')

figure
plot(data.t,data.N,'.')
title('North (m)')

figure
plot(data.t,data.psi,'.')
title('Heading (rad)')

figure
plot(data.t,data.gyro,'.')
title('Gyro (rad/s)')

figure
plot(data.t,data.radar,'.')
title('Radar (m/s)')