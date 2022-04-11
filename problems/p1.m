%% Optimal Estimation - Homework 3 - Problem 1

clear
clc
close all

%% Part A

% Simulation Initialization
t_end = 100;
dt = 0.1;
t = dt:dt:t_end;
numSamps = length(t);

% Model Initialization
A_CL = [0 1; -1 -1.4];
A_D = expm(A_CL*dt);
C = [1 0];
B_W = [0; 1];

% Preallocation
x = zeros(2, numSamps);
y = zeros(numSamps,1);

for i = 1:numSamps-1
    
    % Continuous Simulation
    x_dot = A_CL * x(:,i) + B_W * (2*randn);
    x(:,i+1) = x(:,i) + x_dot * dt;

    % Measurement Simulation
    y(i+1) = x(1,i+1) + randn;
    
end

% Simulation Plotting
figure
plot(t,x(1,:),'.')
hold on
plot(t,y,'.')
title('Position Simulation')
legend('Actual', 'Measurement')
xlabel('Time (s)')
ylabel('Value')

%% Part B
% Q = E[w w']

Q = 4;
Q_D = B_W*Q*B_W'*dt;
R_D = 1;

%% Part C

% Preallocation & Initialization
x_hat = zeros(2, numSamps);
P = eye(2);
priorP = zeros(numSamps,2);
postP = zeros(numSamps,2);

for i = 1:numSamps-1
    
    % Time Update (Propagation - a priori)
    x_hat(:,i+1) = A_D * x_hat(:,i);

    P = A_D * P * A_D' + Q_D;
    priorP(i+1,1) = P(1,1);
    priorP(i+1,2) = P(2,2);

    % Measurement Update (Correction - a posteriori)
    K = [0.1959;0.2094];%P * C' * (C*P*C' + R_D)^-1;

    x_hat(:,i+1) = x_hat(:,i+1) + K * (y(i+1) - C * x_hat(:,i+1));

    P = (eye(2) - K*C) * P;
    postP(i+1,1) = P(1,1);
    postP(i+1,2) = P(2,2);

end

% Steady-State Kalman Gain
K_SS = K;

% Steady-State Covariance
priorSS = priorP(end,:);
postSS = postP(end,:);

% Steady-State Estimator Poles
% KFpoles = eig(K);

% Position Estimate Plotting
figure
plot(t,x(1,:),'.')
hold on
plot(t,x_hat(1,:),'.')
title('Position: Actual vs. Estimate')
legend('Actual','Estimate')
xlabel('Time (s)')
ylabel('Value')

% Covariance Plotting
figure
P1 = plot(t,priorP, 'b');
hold on
P2 = plot(t,postP, 'r');
title('Covariance vs. Time')
legend([P1(1) P2(1)],{'Prior','Posterior'})
xlabel('Time (s)')
ylabel('Covariance (units^2)')

%% Part D

dx = x - x_hat;
std_x = std(dx,0,2);
N1 = vecnorm([std_x(1) std_x(2)]);

%% Part E

Q = 0.79;
Q_D = B_W*Q*B_W'*dt;
R_D = 1;

% Preallocation & Initialization
x_hat = zeros(2, numSamps);
P = eye(2);
priorP = zeros(numSamps,2);
postP = zeros(numSamps,2);

for i = 1:numSamps-1
    
    % Time Update (Propagation - a priori)
    x_hat(:,i+1) = A_D * x_hat(:,i);

    P = A_D * P * A_D' + Q_D;
    priorP(i+1,1) = P(1,1);
    priorP(i+1,2) = P(2,2);

    % Measurement Update (Correction - a posteriori)
    K = P * C' * (C*P*C' + R_D)^-1;

    x_hat(:,i+1) = x_hat(:,i+1) + K * (y(i+1) - C * x_hat(:,i+1));

    P = (eye(2) - K*C) * P;
    postP(i+1,1) = P(1,1);
    postP(i+1,2) = P(2,2);

end

% Steady-State Kalman Gain
K_SS = K;

% Steady-State Covariance
priorSS = priorP(end,:);
postSS = postP(end,:);

% Steady-State Estimator Poles
% KFpoles = eig(K);

% Position Estimate Plotting
figure
plot(t,x(1,:),'.')
hold on
plot(t,x_hat(1,:),'.')
title('Position: Actual vs. Estimate')
legend('Actual','Estimate')
xlabel('Time (s)')
ylabel('Value')

% Covariance Plotting
figure
P1 = plot(t,priorP, 'b');
hold on
P2 = plot(t,postP, 'r');
title('Covariance vs. Time')
legend([P1(1) P2(1)],{'Prior','Posterior'})
xlabel('Time (s)')
ylabel('Covariance (units^2)')

dx = x - x_hat;
std_x = std(dx,0,2);
N2 = vecnorm([std_x(1) std_x(2)]);