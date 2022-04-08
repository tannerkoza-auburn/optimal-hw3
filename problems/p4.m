%% Optimal Estimation - Homework 3 - Problem 4

clear
clc
close all

%% Part A

A = [-2.62 12;
    -0.96 -2];
B = [14; 1];

w = 0.1;
Bw = [1;1];

H = [1 0];

% Simulation Initialization
t_end = 5;
dt = 0.1;
t = dt:dt:t_end;
numSamps = length(t);

X = zeros(2, numSamps);
y = zeros(numSamps, 1);
del = 1;

sys = ss(A,B,H,0);
sysd = c2d(sys,dt);

% Dynamic Simulation
for i = 1:numSamps-1
    xdot = A * X(:,i) + B * del;
    X(:,i+1) = X(:,i) + xdot * dt;
    y(i+1) = H * X(:,i+1) + (0.1*randn);
    x = X(:,i+1);
end

A_D = expm(A*dt);

Q_D = [1 0; 0 1];
R_D = 0.01;


% Preallocation & Initialization
x_hat = zeros(2, numSamps);
P = eye(2);
priorP = zeros(numSamps,2);
postP = zeros(numSamps,2);

for i = 1:numSamps-1
    
    % Time Update (Propagation - a priori)
    x_hat(:,i+1) = A_D * x_hat(:,i) + sys.B * del;

    P = A_D * P * A_D' + Q_D;
    priorP(i+1,1) = P(1,1);
    priorP(i+1,2) = P(2,2);

    % Measurement Update (Correction - a posteriori)
    K = P * H' * (H*P*H' + R_D)^-1;

    x_hat(:,i+1) = x_hat(:,i+1) + K * (y(i+1) - H * x_hat(:,i+1));

    P = (eye(2) - K*H) * P;
    postP(i+1,1) = P(1,1);
    postP(i+1,2) = P(2,2);

end

KFeig = eig(A-K*H);

figure
subplot(2,1,1)
plot(t,X(1,:))
hold on
plot(t,x_hat(1,:))
title('Yaw Rate: Actual & Estimate')
legend('Actual','Estimate')
xlabel('Time (s)')
ylabel('Yaw Rate (rad/s)')

subplot(2,1,2)
plot(t,X(2,:))
hold on
plot(t,x_hat(2,:))
title('Side Slip: Actual & Estimate')
legend('Actual','Estimate')
xlabel('Time (s)')
ylabel('Side Slip (rad)')

clearvars

%% Part B

A = [-2.42 4;
    -0.99 -2];
B = [18; 1];

H = [1 0];

% Simulation Initialization
t_end = 5;
dt = 0.1;
t = dt:dt:t_end;
numSamps = length(t);

X = zeros(2, numSamps);
y = zeros(numSamps, 1);
del = 1;

% Dynamic Simulation
for i = 1:numSamps-1
    xdot = A * X(:,i) + B * del;
    X(:,i+1) = X(:,i) + xdot * dt;
    y(i+1) = H * X(:,i+1) + (0.1*randn);
    x = X(:,i+1);
end

A_A = [-2.62 12;
    -0.96 -2];

A_D = expm(A_A*dt);
Q_D = 5*eye(2);
R_D = 0.01;

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
    K = P * H' * (H*P*H' + R_D)^-1;

    x_hat(:,i+1) = x_hat(:,i+1) + K * (y(i+1) - H * x_hat(:,i+1));

    P = (eye(2) - K*H) * P;
    postP(i+1,1) = P(1,1);
    postP(i+1,2) = P(2,2);

end

KFeig = eig(A-K*H);

figure
subplot(2,1,1)
plot(t,X(1,:))
hold on
plot(t,x_hat(1,:))
title('Yaw Rate: Actual & Estimate')
legend('Actual','Estimate')
xlabel('Time (s)')
ylabel('Yaw Rate (rad/s)')

subplot(2,1,2)
plot(t,X(2,:))
hold on
plot(t,x_hat(2,:))
title('Side Slip: Actual & Estimate')
legend('Actual','Estimate')
xlabel('Time (s)')
ylabel('Side Slip (rad)')

clearvars

%% Part C

A = [-2.42 4;
    -0.99 -2];
B = [18; 1];

H = eye(2);

% Simulation Initialization
t_end = 5;
dt = 0.1;
t = dt:dt:t_end;
numSamps = length(t);

X = zeros(2, numSamps);
y = zeros(2, numSamps);
del = 1;

% Dynamic Simulation
for i = 1:numSamps-1
    xdot = A * X(:,i) + B * del;
    X(:,i+1) = X(:,i) + xdot * dt;
    y(:,i+1) = H * X(:,i+1) + [(0.1*randn);(0.5*randn)];
    x = X(:,i+1);
end

A_A = [-2.62 12;
    -0.96 -2];

A_D = expm(A_A*dt);
Q_D = [5 0;0 0.1];
R_D = [0.01 0;0 0.25];

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
    K = P * H' * (H*P*H' + R_D)^-1;

    x_hat(:,i+1) = x_hat(:,i+1) + K * (y(:,i+1) - H * x_hat(:,i+1));

    P = (eye(2) - K*H) * P;
    postP(i+1,1) = P(1,1);
    postP(i+1,2) = P(2,2);

end

KFeig = eig(A-K*H);

figure
subplot(2,1,1)
plot(t,X(1,:))
hold on
plot(t,x_hat(1,:))
title('Yaw Rate: Actual & Estimate')
legend('Actual','Estimate')
xlabel('Time (s)')
ylabel('Yaw Rate (rad/s)')

subplot(2,1,2)
plot(t,X(2,:))
hold on
plot(t,x_hat(2,:))
title('Side Slip: Actual & Estimate')
legend('Actual','Estimate')
xlabel('Time (s)')
ylabel('Side Slip (rad)')
