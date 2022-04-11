%% Optimal Estimation - Homework 3 - Problem 3

clear 
clc
close all

%% Part A

data = readtable('hw3_3.txt');
data.Properties.VariableNames = {'t','E','N','psi','gyro','v'};

% Time Initialization
numSamps = length(data.t);

% EKF Initialization
xhat = zeros(7,numSamps);
xhat(1,1) = data.E(1);
xhat(2,1) = data.N(1);
xhat(3,1) = data.psi(1);
xhat(4,1) = data.gyro(1);
xhat(5,1) = data.v(1);

procNoise = 0.01;
measNoise = .01;

P = eye(7); 

Q = procNoise * eye(7);
Q(6,6) = 0.00001;
Q(7,7) = 0.00001;

R = measNoise^2 * eye(5);

H = [eye(5) zeros(5,2)];
H(4,6) = 1;
H(5,7) = 1;

P_log = cell(numSamps,1);
P_log{1} = P;

% EKF
for i = 1:numSamps-1

    % Time Step Calculation
    dt = data.t(i+1) - data.t(i);

    % Time Update
    xhat(1,i+1) = xhat(1,i) + xhat(5,i)*sin(xhat(3,i))*dt;
    xhat(2,i+1) = xhat(2,i) + xhat(5,i)*cos(xhat(3,i))*dt;
    xhat(3,i+1) = xhat(3,i) + xhat(4,i)*dt;
    xhat(4,i+1) = xhat(4,i) - xhat(6,i);
    xhat(5,i+1) = xhat(5,i) - xhat(7,i);
    xhat(6,i+1) = xhat(6,i);
    xhat(7,i+1) = xhat(7,i);

    A = eye(7);
    A(1,3) = xhat(5,i)*cos(xhat(3,i))*dt;
    A(1,5) = sin(xhat(3,i))*dt;
    A(2,3) = -xhat(5,i)*sin(xhat(3,i))*dt;
    A(2,5) = cos(xhat(3,i))*dt;
    A(3,4) = dt;
    A(4,6) = -1;
    A(5,7) = -1;

    P = A*P*A' + Q;
    
    % Measurement Update
    
    if i >= 1300
        H = zeros(5,7);
    end
    
    K = P*H'*(H*P*H' + R)^-1;
    
    z = [data.E(i+1); data.N(i+1); data.psi(i+1); data.gyro(i+1);
        data.v(i+1)];

    xhat(:,i+1) = xhat(:,i+1) + K*(z - H*xhat(:,i+1));

    P = (eye(7) - K*H)*P;

    P_log{i+1} = P;

end

figure
plot(data.t,data.E,'.')
hold on
plot(data.t,xhat(1,:),'.')
title('East (m)')

figure
plot(data.t,data.N,'.')
hold on
plot(data.t,xhat(2,:),'.')
title('North (m)')

figure
plot(data.t,data.psi,'.')
hold on
plot(data.t,xhat(3,:),'.')
title('Heading (rad)')

figure
plot(data.t,data.gyro,'.')
hold on
plot(data.t,xhat(4,:),'.')
title('Gyro (rad/s)')

figure
plot(data.t,data.v,'.')
hold on
plot(data.t,xhat(5,:),'.')
title('Radar (m/s)')

figure
plot(data.t,xhat(6,:))
hold on
plot(data.t,xhat(7,:))
title('Mesaurement Biases')
legend('Gyro Bias','Radar Bias')

%% Part B
a = 0;
b = 0;
for i = 1:199
    dt = data.t(i+1) - data.t(i);
    a = a + data.gyro(i+1)*dt;
    b = b + data.v(i+1)*dt;
end
