%% Optimal Estimation - Homework 3 - Problem 2

clear 
clc
close all

%% Part A

% Answer: This happend because the filter is assuming there are no dynamics
% and no process noise. Therefore, when the bias changes, the filter reacts
% poorly.

dt = 0.1;

A_C = 0;
A_D = expm(A_C * dt);
Q_D = 0.1;
H = 1;
R_D = 1;

data = readtable('hw3_2.txt');
data.Properties.VariableNames = {'t','y'};
numSamps = length(data.t);

xhat = zeros(numSamps,1);
x = 0;
P = 1;
K = zeros(numSamps,1);

for i = 1:numSamps
    
    x = A_D * x;
    
    P = A_D * P * A_D' + Q_D;
    
    K(i) = P * H' * (H*P*H' + R_D)^-1;
    
    xhat(i) = x + K(i) * (data.y(i) - H * x);
    
    P = (eye(1) - K(i)*H) * P;
    
    x = xhat(i);

end

figure
plot(data.t, data.y)
hold on
plot(data.t, xhat)
title('Bias')
legend('Measurement','Estimate')

figure
plot(data.t,K)
title('Bias Estimate Kalman Gain')

%% Part B

% Answer: The inclusion of process noise allows for the bias to be tracked
% when it changes. However, the certainty of the state estimate decreases
% (P increases) and doesn't approach 0 as closely.

%% Part C

% Answer: They are the same because 
y_filt = filter(sqrt(Q_D),[1 -(1-sqrt(Q_D))],data.y,data.y(1));

figure
plot(data.t,xhat)
hold on
plot(data.t,y_filt)
title('Kalman vs. Low Pass Filter')
legend('Kalman','Low Pass')
